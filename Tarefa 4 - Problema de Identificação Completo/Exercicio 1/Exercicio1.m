clear;
clc;
close all;

data = load('dados_tarefa4.txt');
%% Carregando dados

t = data(:, 1);
u = data(:, 2);
y = data(:, 3);

un = u - u(1);
yn = (y-mean(y(1:50)));

N = length(t);
Ts = t(2) - t(1);

% Plot inicial dos dados

figure(1);
subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do sistema');

subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada');

% A partir do plot dos dados de entrada e saida, percebe-se um
% comportamento de saida oscilatorio dos dados, e uma boa ideia comecar com um
% modelo de 2a ordem

% ------ 
%% Antes de prosseguirmos, analisaremos a autocorrelacao da saida do sinal  
% para ter certeza que o intervalo de amostragem Ts eh o correto.
% Caso Ts seja alterado, entao os dados de identificacao e validacao
% deverao ser ajustados (PROVAVELMENTE VAO SER....)
%

figure(2);
myccf2(y, N, 0, 1,'k');
xlabel('Atraso');
title('Função de Autocorrelação da saída do sistema');

% Apos ver o grafico, observa-se que o primeiro minimo apos o zero ocorre
% em x=33. De acordo com o teorema da secao 12.2.3 o primeiro minimo deveria ocorrer
% entre o tempo 10 e 20 e esta ocorrendo no tempo 33. 

% Como o valor do minimo esta ocorrendo em 33, entao devemos decimar o
% sinal para que ele fique dentro da faixa estipulada. Faremos uma
% decimacao de fator 2, logo escolheremos 1 a cada 2 amostras iniciais.


%% Fazendo decimacao dos dados de fator 2
% Serao escolhidas 1 a cada 2 amostras

t_decimado_2 = t(1:2:end);
u_decimado_2 = u(1:2:end);
y_decimado_2 = y(1:2:end);

N_dec = length(t_decimado_2); % Novo n decimado
Ts_dec = t_decimado_2(2) - t_decimado_2(1); % Nono Ts decimado

% Agora vamos dividir os dados em 

figure(3);
myccf2(y_decimado_2, N_dec, 0, 1,'k');
xlabel('Atraso');
title('FAC da saída do sistema - Decimada de Fator 2 ');


% Avaliando a figura 3, vemos que a decimacao funcional, pois o novo minimo
% esta em x = 16, ou seja menor que 20.

%% Agora faremos a decimacao dos dados do arquivo, encontrando um novo t, u e y
%  e encontraremos um novo Ts para os dados do sistema.


% Plot dos dados decimados 
figure(4);
subplot(211)
plot(t_decimado_2, y_decimado_2, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do sistema - decimado fator 2');

subplot(212)
plot(t_decimado_2, u_decimado_2);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada - decimado fator 2');


% Em um primeiro momento, serao usados os dados a partir do instante 80 para identificacao
% e do instante 0 a 80 para validacao

% PRE-PROCESSAMENTO - ITEM i
% Dados de identificacao:

t_idt = t_decimado_2(81/Ts_dec:end , 1); 
u_idt = u_decimado_2(81/Ts_dec:end , 1);
y_idt = y_decimado_2(81/Ts_dec:end , 1);

% Frequencia do sistema = 1/Ts_dec = 10Hz


% Dados de validacao:
t_val = t_decimado_2(1:80/Ts_dec , 1); 
u_val = u_decimado_2(1:80/Ts_dec , 1);
y_val = y_decimado_2(1:80/Ts_dec , 1);

%% Identificar se existe atraso puro de tempo nos dados de idetificacao
figure(5);
myccf2([y_idt u_idt], length(y_idt), 1, 1,'k');
xlabel('Atraso');
title('FCC entre o sinal de entrada e saída dados de identificação');

%% Avaliando qual ordem do criterio de akaike usar a partir dos dados de identificacao

data = iddata(y_idt, u_idt, Ts_dec);
akaike = zeros(10,1);
bayes = zeros(10,1);

for i=1:10
    
    a = arx(data, [i i-1 0]);
    akaike(i) = aic(a);
    bayes(i) = aic(a, 'BIC');    
end

figure(6);
plot(akaike);
xlabel('Ordem do modelo')
ylabel('AIC (n\theta)');
title('Critério de informação de Akaike');

figure(7);
plot(bayes);
xlabel('Ordem do modelo')
ylabel('BIC (n\theta)');
title('Critério de informação de Bayes');


% Apos executar o criterio de akaike e Bayes para modelos ARX de 1 a 10 ordem,
% nota-se graficamente, que utilizar um modelo de 3a ordem e a melhor
% escolha de representacao.



%% Criacao do modelo ARX de 3a ordem
% serao usados os dados de identificacao: t_idt, u_idt, y_idt

% Atraso puro de tempo d = 20
d = 20;

% Retirando offset do y(0)


Psi = [y_idt((3+d):(end - 1)), y_idt((2+d):(end-2)), y_idt((d+1):(end-3)), u_idt(1:(end-3-d)), u_idt(2:(end-2-d)), u_idt(3:(end-1-d))];
params = pinv(Psi)*y_idt((4+d): end);

H = tf([params(6) params(5) params(4)], [1 -params(1) -params(2) -params(3)], Ts_dec, 'iodelay', d);

% ------------
y_modelo = lsim(H, u_idt);


%% Validacao do modelo encontrado um Passo a frente

y_estimado = lsim(H, u_val); % Dados de validacao

y_passo_frente = zeros(length(y_val)-3,1);

for i=4:length(y_val)
    y_passo_frente(i-3) = [y_estimado(i-1) y_estimado(i-2) y_val(i-3) u_val(i-3) u_val(i-2) u_val(i-1)]*params;
end

% simulacao um passo a frente
figure(8)
plot(t_val(1:end-3), y_val(1:end-3), 'b');
hold on
plot(t_val(1:end-3), y_passo_frente, 'r--');
title('Validação um passo a frente');
xlabel('Tempo (s)');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');



% indice RMSE
erro_passo_frente = y_val(1:end-3)  - y_passo_frente;
rmse_passo_frente = sqrt(sum(erro_passo_frente.^2)) / sqrt(sum((y_val(1:end-3) -mean(y_val(1:end-3) )).^2));
disp('RMSE simulacao 1 passo a frente');
disp(rmse_passo_frente);


% Residuos do modelo
figure(11);
myccf2(erro_passo_frente, length(erro_passo_frente), 0, 1,'k');
xlabel('Atraso');
title('FAC erro de um passo a frente');


figure(12);
myccf2([u_val(1:end-3) erro_passo_frente], length(erro_passo_frente), 0, 1,'k');
xlabel('Atraso');
title('FCC erro de um passo a frente');



%% Validacao do modelo encontrado usando Simulação livre

yPassoSimulacaoLivre = validacao_livre_3ordem(u_val, y_estimado, params);

figure(13)
plot(t_val(1:end-3), y_val(1:end-3), 'b');
hold on
plot(t_val(1:end-3), yPassoSimulacaoLivre, 'r--');
title('Validação simulação livre');
xlabel('Tempo (s)');
leg = legend({'y','$\hat{y}$'});
set(leg, 'Interpreter', 'latex');


% indice RMSE
erro_simulacao_livre = y_val(1:end-3) - yPassoSimulacaoLivre;

rmse_simulacao_livre = sqrt(sum(erro_simulacao_livre.^2)) / sqrt(sum((y_val(1:end-3)-mean(y_val(1:end-3))).^2));
disp('RMSE simulacao livre');
disp(rmse_simulacao_livre);


% Residuos do modelo
figure(14);
n_val = length(t_val);
myccf2(erro_simulacao_livre, length(erro_simulacao_livre), 0, 1,'k');
xlabel('Atraso');
title('FAC erro de simulacao livre');

figure(15);
myccf2([u_val(1:end-3) erro_simulacao_livre], length(erro_simulacao_livre), 0, 1,'k');
xlabel('Atraso');
title('FCC erro de simulação livre');
