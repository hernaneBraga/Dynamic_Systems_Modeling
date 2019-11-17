clear;
close all;
clc;


%% Letra a) Escolhendo FT tempo continuo e discreto

% Definindo a funcao de transferencia no tempo continuo
qsi = 0.2;
k_eq = 1;
tau = 20;

num = k_eq;
den = [1*tau^2 (2*qsi*tau) 1];
Hs = tf(num,den);

% Convertendo para o tempo discreto
Ts = 8;
Hz = c2d(Hs, Ts);

figure(1)
step(Hs, 'b')
hold on
step(Hz, 'r')
title('Função de trasferência: tempo contínuo (H_s) e discreto (H_z)');
leg = legend({'Hs','Hz'});
set(leg, 'Interpreter', 'latex');

% Confirmando estabilidade dos sistemas:
figure(2)
pzmap(Hz, 'r')
title('Diagrama de polos e zeros da função discretizada');

% Definindo os parametros a1, a2, b1, b2:
[b,a] = tfdata(Hz,'v');
a1 = a(2);
a2 = a(3);
b1 = b(2);
b2 = b(3);


%% Letra a) Gerando sinal PRBS
% Tb do sinal PRBS deve ser:  tau_min/10 <= Tb <= tau_min/3
% De acordo com a secao 4.3.1 do Livro: Introducao a Identificacao de
% sistemas 3a ed.

Tb = 2;
b = 10;
N = 800;
u_prbs = prbs(N, b, Tb);
figure(3);
myccf2(u_prbs', N*2, 1, 1,'k');
xlabel('Atraso');
title('Função de Autocorrelação da entrada PRBS');



%% Letra a) Definindo entradas: Sem Ruido - y1
tam = N;
u1 = u_prbs'; 
t = 0:Ts:(N-1);

y1 = zeros(tam,1);  
for k=3:tam
   y1(k) = -a1*y1(k-1) + (-a2)*y1(k-2) + b1*u1(k-1,1) + b2*u1(k-2,1);
end

%% Letra a) Definindo entradas: Com Ruido - ym

S_R = 25; % Relação Sinal/Ruído em db

e = vetor_erro(y1, S_R); % Retorna o vetor de erro

ym = zeros(tam,1); 
for k=3:tam
   ym(k) = -a1*y1(k-1) + (-a2)*y1(k-2) + b1*u1(k-1) + b2*u1(k-2) + e(k);
end

% % Avaliar se a amostragem do sinal esta ok:
% figure(3333);
% myccf2(y1, N, 0, 1,'k');
% xlabel('Atraso');
% title('Função de Autocorrelação da saída do sistema');
% figure(4444);
% myccf2(ym, N, 0, 1,'k');
% xlabel('Atraso');
% title('Função de Autocorrelação da saída do sistema');



% Comparacao da saída real e da saída com ruído S/R [db]
figure(4);
subplot(211)
plot(y1, 'b') 
hold on
plot(ym, 'r.-.') 
title(['Comparação de saída sem ruído e relação S/R =' num2str(S_R) 'db']);
legend('y',['y_m' num2str(S_R) 'db']);


subplot(212)
plot(u_prbs);
xlabel('Amostras')
ylabel('Amplitude');
title('Entrada PRBS');

%% Letra c) I - Modelo ARX e simulacoes sem ruido - Estimacao do modelo

% Definindo o modelo ARX de segunda ordem sem ruido
Psi1=[y1(2:N-1) y1(1:N-2) u1(2:N-1) u1(1:N-2)];
params = pinv(Psi1)*y1(3:end);

% Validando a resposta encontrada
Hz_limpo = tf([params(3) params(4)], [1 -params(1) -params(2)], Ts);

disp('Parametros encontrados sem ruído:');
disp(params);

% Diagrama de polos e zeros
figure(5);
pzmap(Hz, 'b')
hold on
pzmap(Hz_limpo, 'r')
title('Diagrama de polos e zeros');
leg = legend({'Hz','$\hat{Hz}$'});
set(leg, 'Interpreter', 'latex');

% Pela figura, percebe-se que os zeros e polos estao no mesmo lugar

%% Letra c) I - Modelo ARX e simulacoes sem ruido - Validaco passo frente e livre

u2 = ones(1, N);
yreal = step(Hz,t);

% =============== SIMULACAO PASSO A FRENTE ======================
y_passo_frente = validacao_passo_frente(u2, yreal, 2, params);
y_passo_frente = y_passo_frente';

figure(6);
plot(yreal, 'b') 
hold on
plot(y_passo_frente, 'r.-.') 
title('Validação um passo a frente: Resposta ao degrau');
xlabel('Amostras');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');

rmse_passo_frente = RMSE(yreal, y_passo_frente);
disp('RMSE simulacao de um passo a frente sem ruído:'); 
disp(rmse_passo_frente);


% =============== SIMULACAO LIVRE ======================
y_livre = validacao_livre_2ordem(u2, yreal , params);

figure(7);
plot(yreal, 'b') 
hold on
plot(y_livre, 'r.-.') 
title('Validação simulação livre: Resposta ao degrau');
xlabel('Amostras');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');

rmse_livre = RMSE(yreal(1:end-3), y_livre);
disp('RMSE simulacao livre:'); 
disp(rmse_livre);

% =============== ANÁLISE DOS RESÍDUOS - SEM RUÍDO ======================
% Mostrando a autocorrelação do vetor de residuos das duas validaçoes
figure(71)
subplot(211)
myccf2((yreal - y_passo_frente), length(y_passo_frente), 0, 1,'k');
xlabel('Atraso');
title('FAC resíduos de um passo a frente');

subplot(212)
myccf2((yreal(1:end-3) - y_livre), length(y_livre), 0, 1,'k');
xlabel('Atraso');
title('FAC resíduos da simulacao livre');

%% Letra c) II - Modelo ARX e simulacoes com ruido

% Definindo o modelo ARX de segunda ordem a partir da saida com ruido
Psi1_ruido=[ym(2:N-1) ym(1:N-2) u1(2:N-1) u1(1:N-2)];
params_ruido = pinv(Psi1_ruido)*ym(3:end);

% Validando a resposta encontrada
Hz_ruido = tf([params_ruido(3) params_ruido(4)], [1 -params_ruido(1) -params_ruido(2)], Ts);

disp(['Parametros encontrados com Sinal/Ruído = ' num2str(S_R) 'db']);
disp(params_ruido);

% --- Diagrama de polos e zeros
figure(8);
pzmap(Hz, 'b')
hold on
pzmap(Hz_ruido, 'r')
title(['Diagrama de polos e zeros: S/R = ' num2str(S_R) 'db']);
legend('Hz',['Hz ' num2str(S_R) 'db']);

%% Letra c) II - Modelo ARX e simulacoes com ruido - Validaco passo frente e livre

ym_hat = step(Hz_ruido, t); % Saida estimada ao degrau e ruido

% =============== SIMULACAO PASSO A FRENTE - RUÍDO ======================

y_passo_frente_ruido = validacao_passo_frente(u2, ym_hat, 2, params_ruido);
y_passo_frente_ruido = y_passo_frente_ruido';

figure(9);
plot(t,yreal, 'b') 
hold on
plot(t,y_passo_frente_ruido, 'r.-.') 
title(['Validação um passo a frente: S/R = ' num2str(S_R) 'db']);
xlabel('Amostras');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');

rmse_passo_frente = RMSE(yreal, y_passo_frente_ruido);
disp(['RMSE simulacao de um passo a frente: S/R = ' num2str(S_R) 'db']); 
disp(rmse_passo_frente);


% =============== SIMULACAO LIVRE - RUÍDO ======================

y_livre_ruido = validacao_livre_2ordem(u2, ym_hat , params_ruido);

figure(88);
%plot(t(1:end-3), yreal(1:end-3), 'b') 
plot(yreal, 'b') 
hold on
plot( y_livre_ruido, 'r.-.') 
%plot(t(1:end-3), y_livre_ruido, 'r.-.') 
title(['Validação simulação livre: S/R = ' num2str(S_R) 'db']);
xlabel('Amostras');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');

rmse_livre = RMSE(yreal(1:end-3), y_livre_ruido);
disp('RMSE simulacao livre:'); 
disp(rmse_livre);



% =============== ANÁLISE DOS RESÍDUOS ======================

figure(91)
subplot(211)
myccf2((yreal - y_passo_frente_ruido), length(y_passo_frente), 0, 1,'k');
xlabel('Atraso');
title(['FAC resíduos de um passo a frente: S/R = ' num2str(S_R) 'db']);

subplot(212)
myccf2((yreal(1:end-3) - y_livre_ruido), length(y_livre_ruido), 0, 1,'k');
xlabel('Atraso');
title(['FAC resíduos da simulação livre: S/R = ' num2str(S_R) 'db']);

%% Letra d) Avaliando o efeito do aumento da relação sinal/ruído

%Definindo a saida com novos niveis de relação S/R
ruidos = [5 10 15 20];

%Matrizes que guardarão os resultados
params_finais = zeros(length(ruidos) , 4);
y_p_frente = zeros(length(ruidos) , length(y_passo_frente) );
rmse_p_frente = zeros(4, 1);
y_p_livre = zeros(length(ruidos) , length(y_livre_ruido) );
rmse_p_livre = zeros(4 , 1 );

% Calculando: y_passo_frente, y_simulacao_livre, parametros e rmse
for i=1:length(ruidos)
    e = vetor_erro(y1, ruidos(i)); % Gerando erro
    ym = zeros(tam,1); % Gerando saida com ruido de processo
    
    for k=3:tam
        ym(k) = -a1*y1(k-1) + (-a2)*y1(k-2) + b1*u1(k-1) + b2*u1(k-2) + e(k);
    end
    
    Psi1_ruido=[ym(2:N-1) ym(1:N-2) u1(2:N-1) u1(1:N-2)];  % Matriz de regressores
    params_ruido = pinv(Psi1_ruido)*ym(3:end);
    params_finais(i,:) = params_ruido;
    
    % ---- Passo a frente:
    y_p_frente(i,:) = validacao_passo_frente(u2, ym_hat, 2, params_ruido);
    y_p_livre(i,:) = validacao_livre_2ordem(u2, ym_hat , params_ruido);
    
    rmse_passo_frente(i) = RMSE(yreal, y_p_frente(i,:)');
    rmse_p_livre(i) = RMSE(yreal(1:end-3), y_p_livre(i,:)');
    

end


% ============================== Plot dos resultados:
% ======= Passo a frente
figure(999);
plot(yreal, 'b') 
hold on
plot(y_p_frente(1,:), 'r--') 
hold on
plot(y_p_frente(2,:), 'g--')
hold on
plot(y_p_frente(3,:), 'm--')
hold on
plot(y_p_frente(4,:), 'k--')
title('Validação um passo a frente: Comparação entre diferentes S/R');
xlabel('Amostras');
legend('y',['ym ' num2str(ruidos(1)) 'db'],['ym ' num2str(ruidos(2)) 'db'], ['ym ' num2str(ruidos(3)) 'db'], ['ym ' num2str(ruidos(4)) 'db']);

% ============================== Plot dos resultados:
% ======= SImulação livre  

figure(222);
plot(t(1:end-3), yreal(1:end-3), 'b') 
hold on
plot(t(1:end-3), y_p_livre(1,:), 'r--') 
hold on
plot(t(1:end-3), y_p_livre(2,:), 'g--')
hold on
plot(t(1:end-3), y_p_livre(3,:), 'm--')
hold on
plot(t(1:end-3), y_p_livre(4,:), 'k--')
title('Validação simulação livre: Comparação entre diferentes S/R');
xlabel('Amostras');
legend('y',['ym ' num2str(ruidos(1)) 'db'],['ym ' num2str(ruidos(2)) 'db'], ['ym ' num2str(ruidos(3)) 'db'], ['ym ' num2str(ruidos(4)) 'db']);

% ============================== RMSE:
% Passo a frente
disp('RMSE simulacao de um passo a frente para S/R db:'); 
disp(ruidos);
disp(rmse_passo_frente);

% Simulação livre
disp('RMSE simulacao livre para S/R db:'); 
disp(ruidos);
disp(rmse_p_livre');

% ============================== Diagrama de polos e zeros Hz, 5db e 20db:
Hz_ruido_5db = tf([params_finais(1,3) params_finais(1,4)], [1 -params_finais(1,1) -params_finais(1,2)], Ts);
Hz_ruido_20db = tf([params_finais(4,3) params_finais(4,4)], [1 -params_finais(4,1) -params_finais(4,2)], Ts);


figure(103);
pzmap(Hz, 'b')
hold on
pzmap(Hz_ruido_5db, 'r')
hold on
pzmap(Hz_ruido_20db, 'k')
title('Comparação no diagrama de polos e zeros');
legend('Hz',['Hz ' num2str(ruidos(1)) 'db'], ['Hz ' num2str(ruidos(4)) 'db']);


