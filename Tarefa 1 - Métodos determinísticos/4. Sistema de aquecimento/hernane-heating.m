%% Carregando dados amostrados e definindo parametros iniciais
clear all; clc
dados = load('heating_system.dat');

% 1 parte
t = round(dados(1:376,1)); t0 = t(1);
u = dados(1:376,2); u0 = u(1);
y = dados(1:376,3); y0 = y(1);
u = u-6;

% tudo
% t = round(dados(:,1)); t0 = t(1);
% u = dados(:,2); u0 = u(1);
% y = dados(:,3); y0 = y(1);


y_ajustado = y - y0;
tam = length(t);

%% Plot dos dados carregados

figure(1);
subplot(211); plot(t,y); xlabel('Tempo (s)'); ylabel('y(t) [Cº]'); title('Sistema de calor');
subplot(212); plot(t, u,'k-.','LineWidth',2); xlabel('Tempo (s)');ylabel('u(t)]V]');

%% Cálculo do ganho K
delta_u = (u(end) - u(1));
delta_y = y_ajustado(end)- y_ajustado(1);
K = delta_y/delta_u;
K = 2.75*K;
%K = K*2;
%% Resposta complementar:

[~, t35] = min(abs(y_ajustado - y_ajustado(end)*0.35));
[~, t85] = min(abs(y_ajustado - y_ajustado(end)*0.85));

tau1 = 0.682*(t(t85)-t(t35));
theta1 = 0;


G1a = tf(K, [tau1 1], 'ioDelay', theta1);
resp = lsim(G1a, u, t) + y0;

figure(3);
plot(t,y); xlabel('Tempo (s)'); ylabel('y(t) [Cº]'); title('Função de Transferência');
hold on
plot(t,resp,'r-.','LineWidth',2); xlabel('Tempo (s)');ylabel('G(t)');
