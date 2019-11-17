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

%% Letra c) I - Modelo ARX e simulacoes sem ruido - Estimacao do modelo de 1 ordem

% Definindo o modelo ARX de primeira ordem sem ruido
Psi1=[y1(3:N-1) y1(2:N-2) y1(1:N-3) u1(3:N-1) u1(2:N-2) u1(1:N-3)];
params = pinv(Psi1)*y1(4:end);

Hz_limpo = tf([params(4) params(5) params(6)], [1 -params(1) -params(2) -params(3)], Ts);

disp('Parametros encontrados sem ruído:');
disp(params);

%% Letra c) I - Modelo ARX e simulacoes sem ruido - Validaco passo frente e livre

u2 = ones(1, N);
yreal = step(Hz,t);
y3ordem = step(Hz_limpo,t);

% =============== SIMULACAO PASSO A FRENTE ======================
y_passo_frente = validacao_passo_frente(u2, yreal, 3, params);
y_passo_frente = y_passo_frente';

rmse_passo_frente = RMSE(yreal, y_passo_frente);
disp('RMSE simulacao de 3º ordem um passo a frente sem ruído:'); 
disp(rmse_passo_frente);


% ===============  SIMULACAO LIVRE ======================

y_livre = validacao_livre_3ordem(u2, yreal, params);

rmse_livre = RMSE(yreal(1:end-3), y_livre);
disp('RMSE simulacao livre sem ruído (3a ordem):'); 
disp(rmse_livre);


% =============== ANÁLISE DOS RESÍDUOS - SEM RUÍDO ======================
% Mostrando a autocorrelação do vetor de residuos das duas validaçoes
figure(71)
subplot(211)
myccf2((yreal - y_passo_frente), length(y_passo_frente), 0, 1,'k');
xlabel('Atraso');
title('FAC resíduos de um passo a frente - 3º ordem (sem ruído)');

subplot(212)
myccf2((yreal(4:end) - y_livre), length(y_livre), 0, 1,'k');
xlabel('Atraso');
title('FAC resíduos da simulacao livre - 3º ordem (sem ruído)');

%% Letra c) II - Modelo ARX e simulacoes com ruido

% Resposta real definida pela saida medida de 2a ordem:
Psi1_ruido_real=[ym(2:N-1) ym(1:N-2) u1(2:N-1) u1(1:N-2)];
params_ruido_real = pinv(Psi1_ruido_real)*ym(3:end);
Hz_ruido_real = tf([params_ruido_real(3) params_ruido_real(4)], [1 -params_ruido_real(1) -params_ruido_real(2)], Ts);

% Definindo o modelo ARX de segunda ordem a partir da saida com ruido
Psi1_ruido=[ym(3:N-1) ym(2:N-2) ym(1:N-3) u1(3:N-1) u1(2:N-2) u1(1:N-3)];
params_ruido = pinv(Psi1_ruido)*ym(4:end);
Hz_ruido = tf([params_ruido(4) params_ruido(5) params_ruido(6)], [1 -params_ruido(1) -params_ruido(2) -params_ruido(3)], Ts);

disp(['Parametros encontrados de 3º ordem com Sinal/Ruído = ' num2str(S_R) 'db']);
disp(params_ruido);

%% Letra c) II - Modelo ARX e simulacoes com ruido - Validaco passo frente e livre

ym_hat = step(Hz_ruido, t); % Saida estimada ao degrau e ruido
ym_real = step(Hz_ruido_real, t); 

% =============== SIMULACAO PASSO A FRENTE - RUÍDO ======================
y_passo_frente_ruido = validacao_passo_frente(u2, ym_real, 3, params_ruido);
y_passo_frente_ruido = y_passo_frente_ruido';

rmse_passo_frente = RMSE(yreal, y_passo_frente_ruido);
disp(['RMSE simulacao 3º ordem de um passo a frente: S/R = ' num2str(S_R) 'db']); 
disp(rmse_passo_frente);

% =============== FAZER FUNCIONAR - SIMULACAO LIVRE ======================

y_livre_ruido = validacao_livre_3ordem(u2, ym_real, params_ruido);

rmse_livre_uido = RMSE(yreal(1:end-3), y_livre_ruido);
disp(['RMSE simulacao livre: S/R = ' num2str(S_R) 'db']); 
disp(rmse_livre_uido);


% =============== ANÁLISE DOS RESÍDUOS ======================

figure(91)
subplot(211)
myccf2((yreal - y_passo_frente_ruido), length(y_passo_frente), 0, 1,'k');
xlabel('Atraso');
title(['FAC de 3º ordem resíduos de um passo a frente: S/R = ' num2str(S_R) 'db']);

subplot(212)
myccf2((yreal(4:end) - y_livre_ruido), length(y_livre_ruido), 0, 1,'k');
xlabel('Atraso');
title(['FAC de 3º ordem resíduos da simulação livre: S/R = ' num2str(S_R) 'db']);

%% Resultados 

% SIMULAÇÃO UM PASSO A FRENTE
figure(222);
plot(yreal, 'b') 
hold on
plot(y_passo_frente, 'k--') 
hold on
plot(y_passo_frente_ruido, 'r--') 
title('Validação 3º ordem um passo a frente: Resposta ao degrau');
xlabel('Amostras');
leg = legend({'y','$\hat{y}$ sem ruido','$\hat{y}$ 25db'});
set(leg, 'Interpreter', 'latex');


% SIMULAÇÃO LIVRE

figure(333);
plot(yreal, 'b') 
hold on
plot(y_livre, 'k--') 
hold on
plot(y_livre_ruido, 'r--') 
title('Validação 3º ordem simulação livre: Resposta ao degrau');
xlabel('Amostras');
leg = legend({'y','$\hat{y}$ sem ruido','$\hat{y}$ 25db'});
set(leg, 'Interpreter', 'latex');



% Diagrama de polos e zeros
figure(888);
pzmap(Hz, 'b')
hold on
pzmap(Hz_limpo, 'k')
hold on
pzmap(Hz_ruido, 'r')
title('Diagrama de polos e zeros: 2º e 3º ordem');
legend('Hz', 'Hz sem ruído',['Hz ' num2str(S_R) 'db']);
