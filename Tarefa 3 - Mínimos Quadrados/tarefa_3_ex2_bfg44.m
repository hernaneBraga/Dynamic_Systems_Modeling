clear;
clc;
close all;

data = load('BFG44.DAT');
%% Carregando dados

t = data(:, 1);
u = data(:, 2);
y = data(:, 3);

un = u - u(1);
yn = (y-mean(y(1:50)));

N = length(t);
Ts = t(2) - t(1);

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


figure(2);
subplot(211)
plot(t, yn, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta do sistema - normalizada');

subplot(212)
plot(t, un);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada');



%% Hernane
d = 4;
Psi = [y(1+d:N-1), u(1:N-1-d)];
theta = (Psi' * Psi) \ Psi' * y(2+d:N);

tau_hat = - Ts / (theta(1) - 1);
K_hat = (tau_hat * theta(2)) / Ts;

txt = 'Constante de tempo estimada: ';
disp(txt);
disp(tau_hat);
txt = 'Ganho estimado: ';
disp(txt);
disp(K_hat);

%K_hat = 0.35;

G_modelo = tf([K_hat], [tau_hat 1],'ioDelay', (d-1));
y_modelo = lsim(G_modelo, un, t);


H1 = tf(1.338, [3.586 4.459 1], 'ioDelay', 1.9);
H2 = tf(0.0182, [1 0.1824 0.052], 'ioDelay', 4.7);
y_h1 = lsim(c2d(H1, 0.1), un);

figure(3);
plot(t, yn, 'b');
hold on
plot(t, y_h1, 'g');
hold on
plot(t, y_modelo, 'r');
legend('y(s)','H_2(s)','Modelo');
title('y(s) x H_2(s) x Modelo MQ');
xlabel('Tempo (s)');
ylabel('Amplitude');


figure(4);
step(H2, 'b');
hold on
step(G_modelo, 'r-')
legend('H_2(s)','Modelo');
title('H2(s) x Modelo encontrado - Resposta ao impulso');
xlabel('Tempo (s)');
ylabel('Amplitude');


