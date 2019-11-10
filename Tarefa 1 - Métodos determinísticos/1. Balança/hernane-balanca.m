%% Carregando dados amostrados e definindo 
clear all; clc
dados_prato_1 = load('DS100g_prato_1.txt');
dados_prato_2 = load('DS100g_prato_2.txt');
degrau_inicial = 100;

%% Plot dos dados amostrados da balança 1
figure(1);plot(dados_prato_1(:,1),dados_prato_1(:,2));
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao degrau - Prato 1');

%% Tratamento de dados das amostras
posicao_inicial = find(dados_prato_1(:,1) == 0);
posicao_inicial_2 = find(dados_prato_2(:,1) == 0);
y0_dados = dados_prato_1(posicao_inicial,2); 

%Definindo o vetor da amplitude de saida e retidando o valor da amplitude inicial y0
%Assim amplitude inicial do sinal será: y(0) = 0
balanca_1_ajustado = dados_prato_1(posicao_inicial:end, 2) - y0_dados; 

%definicao do vetor de tempo t
t = dados_prato_1(posicao_inicial:end, 1); 
t2 = dados_prato_2(posicao_inicial:end, 1);
tam = length(t);

%% Plot dos dados da balança 1 ajustados iniciando em t = 0 e y(0) = 0
figure(2);plot(t, balanca_1_ajustado);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao degrau ajustada - Prato 1');

%% Definição dos parâmetros da Eq. de Transferência à partir dos dados ajustados: K, Wn e Zeta
% Valor de K:
  %Considerando que um peso de 100g foi colocado no instante zero, e que no final do sinal o valor tende a 0.6 
  %então o ganho K é: K = 0.6/100 -> K = 0,006
    k = balanca_1_ajustado(tam)/degrau_inicial;
    
% Valor de Wn:
  %Usaremos um modelo de 2º ordem, assim Wn pode ser estimado através da fórmula: períodos 
  %Wn = (2pi/T) onde T é o valor de um periodo. 
  %Usando a média de 5 períodos (colhidos visualmente) T~= 0.05
    wn = (2*pi)/0.05;
 
 % Valor de Zeta:
  %Usando a fórmula do material de referência zeta = 0.6/n_ciclos
  %Onde contando visualmente o número de ciclos é aproximadamente 60
    zeta = 0.6/60;
  
%% Definição da função de transferência H(s)
numerador = [k*(wn)^2];

d1 = 1;
d2 = 2*zeta*wn;
d3 = wn^2;
denominador = [d1 d2 d3];

func_trasnf = tf(numerador, denominador);

%% Plot da resposta ao degrau dos dados da balança 1 x função de transferência

u0 = 0;
u = [u0; degrau_inicial*ones(tam - 1 ,1)];
resp = lsim(func_trasnf, u, t) + y0_dados; 

figure(3);plot(t,resp,'r');
hold on
plot(t, dados_prato_1(posicao_inicial:end, 2));
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do prato 1 x Função de transferência ');

%% Plot da resposta ao degrau dos dados da balança 2 x função de transferência
figure(4);plot(t,resp,'r');
hold on
plot(t2, dados_prato_2(posicao_inicial_2:end, 2));
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do prato 2 x Função de transferência ');


