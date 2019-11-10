%% Carregando dados amostrados e definindo parametros iniciais
clear all; clc;
dados = load('torneira3.txt');

u = dados(:,2); u0 = u(1);
y = dados(:,1); y0 = y(1);

y_ajustado = y - y0;


t_amostra = 1; %tempo de amostragem (s) de acordo com o readme.txt
t = [0:t_amostra:(length(y)-1)]';
tam = length(t);

%% Plot dos dados carregados

figure(1);
subplot(211); plot(t,y); xlabel('Tempo (s)'); ylabel('y(t)'); title('Torneira 3');
subplot(212); plot(t, u,'k-.','LineWidth',2); xlabel('Tempo (s)');ylabel('u(t)');

%% Cálculo do ganho K
delta_u = abs(u(end) - u(1));
delta_y = abs(mean(y_ajustado(round(tam*0.66):end))- mean(y_ajustado(1:round(tam*0.3))));
K = delta_y/delta_u;

%% Resposta complementar:

t_10 = round(0.1*tam);
t_33 = round(0.15*tam);


yy = log(abs(1 - y_ajustado./(K * u)));
%figure(2); plot(t, yy); xlabel('t (s)'); ylabel('ln(y(t)/(K*u(t)))');
coef = polyfit(t(t_10:t_33),yy(t_10:t_33),1);
%hold on; plot(t, coef(1)*t + coef(2), 'm-.', 'LineWidth', 2);
tau2a = -1/coef(1);

%%%%%%%
t_10_2 = round(0.1*t_10);
t_33_2 = round(((t_33-t_10)/100)*tam);

yy2 = log(abs(exp(coef(2))*exp(-(t)./tau2a) - (1 - y_ajustado./(K*u))));
coef2 = polyfit(t(t_10_2:t_33_2),yy2(t_10_2:t_33_2),1);
%figure(3);  plot(t, yy2); xlabel('t (s)'); ylabel('v(t)');
%hold on; plot(t, coef2(1)*t + coef2(2), 'm-.', 'LineWidth', 2);
tau2b = -1/coef2(1);

% Atraso de tempo é de aproximadamente 20s (inspecionando visualmente)
theta = 15;

K = 1.85*K;
G2a = tf(K, [tau2a*tau2b  tau2a+tau2b  1], 'ioDelay', theta);
resp = lsim(G2a, u, t) + y0;

% figure(3);
% plot(t,y); xlabel('Tempo (s)'); ylabel('y(t)'); title('Torneira 3');
% hold on
% plot(t,resp,'r-.','LineWidth',2); xlabel('Tempo (s)');ylabel('G(t)');


figure(3);
subplot(211); plot(t,y); xlabel('Tempo (s)'); ylabel('y(t)'); title('Torneira 3');
hold on
plot(t,resp,'r-.','LineWidth',2); xlabel('Tempo (s)');ylabel('G(t)');
subplot(212); plot(t, u,'k-.','LineWidth',2); xlabel('Tempo (s)');ylabel('u(t)');


%% Testando a FT para a torneira 4
dados = load('torneira4.txt');

u = dados(:,2); u0 = u(1);
y = dados(:,1); y0 = y(1);

y_ajustado = y - y0;

t_amostra = 1; %tempo de amostragem (s) de acordo com o readme.txt
t = [0:t_amostra:(length(y)-1)]';
tam = length(t);

theta = 10;
G2a = tf(K, [tau2a*tau2b  tau2a+tau2b  1], 'ioDelay', theta);
resp = lsim(G2a, u, t) + y0;

figure(4);
subplot(211); plot(t,y); xlabel('Tempo (s)'); ylabel('y(t)'); title('Torneira 4');
hold on
plot(t,resp,'r-.','LineWidth',2); xlabel('Tempo (s)');ylabel('G(t)');
subplot(212); plot(t, u,'k-.','LineWidth',2); xlabel('Tempo (s)');ylabel('u(t)');

%% Testando a FT para a torneira 5

dados = load('torneira5.txt');

u = dados(:,2); u0 = u(1);
y = dados(:,1); y0 = y(1);

y_ajustado = y - y0;

t_amostra = 1; %tempo de amostragem (s) de acordo com o readme.txt
t = [0:t_amostra:(length(y)-1)]';

tam = length(t);

K = K/1.5;
theta = 65; 
G2a = tf(K, [tau2a*tau2b  tau2a+tau2b  1], 'ioDelay', theta);
resp = lsim(G2a, u, t) + y0;

figure(5);
subplot(211); plot(t,y); xlabel('Tempo (s)'); ylabel('y(t)'); title('Torneira 5');
hold on
plot(t,resp,'r-.','LineWidth',2); xlabel('Tempo (s)');ylabel('G(t)');
subplot(212); plot(t, u,'k-.','LineWidth',2); xlabel('Tempo (s)');ylabel('u(t)');