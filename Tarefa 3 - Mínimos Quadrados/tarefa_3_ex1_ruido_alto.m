clc;
clear;
close all;

%% Sistema verdadeiro G(s)
tau = 30; %constante de tempo
K = 5;
theta = 50; %atraso de tempo
Gs = tf([K], [tau 1],'ioDelay', theta);

% Plot da resposta ao impulso
figure(1)
step(Gs, 'b')
title('Resposta ao impulso');

%% Dados de entrada PRBS:
Ts = tau/10;  % Intervalo de amostragem
N = 300;
t = [0:Ts:(N-1)*Ts]';
t2 = [0:Ts:(N-1)*Ts]';
Tb = tau/5; % tempo de permanência em s
m = Tb/Ts; % tempo de permanência em amostras
u = prbs(N, 9,m);


r = randn(N, 1) * 1.5;
u = r'+u;

%%
% Simulação da eq. criada:
yreal = lsim(Gs, u, t);

%Plot da entrada PRBS x saída
figure(2);
subplot(211)
plot(t, yreal, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do sistema - alto ruído');


subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada PRBS - alto ruído');
ylim([-4 4])


%% Estimação do atraso puro de tempo via FCC

figure(3);
[t,ruu,l,B]=myccf2([yreal u'],N,1,1,'k');
xlabel('Atraso')
title('Autocorrelação entre y(k) e u(k)');


%% Metodo dos Minimos Quadrados

d = 17; %Atraso estimado via FCC. 

u_k = u(1:N-d);
y_k = yreal(d+1:N); %y(k)
y_k_1 = yreal(d+1:N); %y(k-1)

% Montando o vetor de atraso y(k-1)
for i = 1:N-d
    if i == 1
        y_k_1(i) = yreal(d);
    else
        y_k_1(i) = y_k(i-1);
    end
end

% Montando a matriz ARX
Psi = [y_k_1, u_k'];
theta = (Psi' * Psi) \ Psi' * y_k;

% Estimando a resposta 
tau_hat = - Ts / (theta(1) - 1); % constante de tempo estimada
K_hat = (tau_hat * theta(2)) / Ts; % ganho estimado

txt = 'Constante de tempo estimada: ';
disp(txt);
disp(tau_hat);
txt = 'Ganho estimado: ';
disp(txt);
disp(K_hat);

%% Montando a nova FT e plotando o gráfico de comparação


G_modelo = tf([K_hat], [tau_hat 1],'ioDelay', (d-1)*Ts);
y_modelo = lsim(G_modelo, u, t2);

figure(4)
subplot(211)
plot(t2, yreal, 'b');
hold on
plot(t2, y_modelo, 'g--');
legend('G(s)','Modelo');
title('G(s) x Modelo encontrado - Entrada PRBS');
xlabel('Tempo (s)');
ylabel('Amplitude');


subplot(212)
step(Gs, 'b');
hold on
step(G_modelo, 'g--')
legend('G(s)','Modelo');
title('G(s) x Modelo encontrado - Resposta ao impulso com ruído alto');
xlabel('Tempo (s)');
ylabel('Amplitude');

