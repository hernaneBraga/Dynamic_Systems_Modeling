%% Definição de variáveis globais e parametros de simulacao
clear all; clc

%% Carregando dados amostrados e definindo 
dados = load('ENS_25.DAT');

% Plot dos dados amostrados da balança 1
figure(1);plot(dados(:,1),dados(:,2));
xlabel('Tempo (s)');
ylabel('Nivel em m');
title('Resposta ao degrau - Bomba 25');

%% Tratamento de dados das amostras
posicao_inicial = find(dados(:,1) == 0);
y0_dados = dados(posicao_inicial,2); 

%Definindo o vetor da amplitude de saida e retidando o valor da amplitude inicial y0
%Assim amplitude inicial do sinal será: y(0) = 0
y_ajustado = dados(posicao_inicial:end, 2) - y0_dados; 

%definicao do vetor de tempo t
t = dados(posicao_inicial:end, 1); 
tam = length(t);
t_amostra = t(2) - t(1); %Intervalo de amostragem

% Plot dos dados da balança 1 ajustados iniciando em t = 0 e y(0) = 0
figure(2);plot(t, y_ajustado);
xlabel('Tempo (s)');
ylabel('Nivel em m');
title('Resposta ao degrau ajustada - Bomba 25');

%% Criação do degrau de entrada
tempo_u = 50; %Duracao do degrau
u = zeros(tam, 1);

%Degrau positivo de 16,34mA a 17,05mA
u(1:tempo_u) = 16.34; 
u(tempo_u:tam) = 17.05;

%% Ganho K
delta_u = u(end)-u(1);
delta_y = y_ajustado(end) - y_ajustado(1);
K = delta_y/delta_u;

%% Normalização de yn e un para o método das áreas
aux = 500;
yn = (y_ajustado - mean(y_ajustado(1:aux))) / (mean( y_ajustado(end-aux:end)) - mean(y_ajustado(1:aux)) );
un = (u - u(1))/(delta_u);

yn = yn(tempo_u+1:end);
un = un(tempo_u+1:end);
tn = t(tempo_u+1:end) - t(tempo_u);


%% Método das áreas
area = sum(t_amostra*(un - yn));

%tau1 = exp(1)*sum(Ts*yn(1:find(t==round(area))));
tau = exp(1)*sum(t_amostra*yn(1:find(round(tn)==round(area))));
theta = area - tau;

% K = 0.2837
%K = K/10; ajuste 1
K = K/20; % ajuste 2
func_trasnf = tf(K, [tau 1], 'ioDelay', theta);

%% Plot da planta 25 x função de transferência
%u0 = 0;
%u = [u0; degrau_inicial*ones(tam - 1 ,1)];

resp = lsim(func_trasnf, u, t) + y0_dados; 
figure(3);plot(t,resp,'r');
hold on
plot(t, dados(posicao_inicial:end, 2), 'b');
xlabel('Tempo (s)');
ylabel('Nivel em m');
title('Resposta da bomba 25 x Função de transferência modelada');

%% Valindo FT com os dados da planta 26

dados = load('ENS_26.DAT');
posicao_inicial = find(dados(:,1) == 0);
y0_dados = dados(posicao_inicial,2); 
t = dados(posicao_inicial:end, 1); 

%Degrau negativo de 17,05mA a 16,34mA
u2(1:tempo_u) = -17.05;
u2(tempo_u:tam) = -16.34;

resp2 = lsim(func_trasnf, u2, t) + y0_dados; 
figure(4);plot(t,resp2,'r');
hold on
plot(t, dados(:, 2), 'b');
xlabel('Tempo (s)');
ylabel('Nivel em m');
title('Resposta da bomba 26 x Função de transferência modelada');


