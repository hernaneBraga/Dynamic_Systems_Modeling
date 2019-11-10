% Hernane Braga Pereira - UFMG Agosto 2019
% Exerc�cio 3.10 - Livro Introdu��o � Identifica��o de Sistemas 3.Ed

% Implementa��o das rotinas usadas s�o de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Dispon�veis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clear
close all

% Definicao da funcao de transferencia
H = tf(1, [1 0.2 0.8]);

% Parametros usados:
N = 15; % N�mero de amostras
t = [0:N-1]'; 

%% 1) Entrada usando impulso
u = zeros(1,N);
u(1)= 1;
y = lsim(H, u, t);

% Estimando saida usando convolucao
diagonal(1:N, 1:N) = 1;
diagonal = diag(diag(diagonal));
y_est = (diagonal^-1)*y;

% Plot dos resultados
figure(1);
subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada: impulso');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');

%% 1.1) Entrada usando impulso - com ruido
u = zeros(1,N);
u(1)= 1;
y = lsim(H, u, t);

% Adicionando ruido a saida
y = y + 0.05*randn(length(y),1);

% Estimando saida usando convolucao
diagonal(1:N, 1:N) = 1;
diagonal = diag(diag(diagonal));
y_est = (diagonal^-1)*y;

% Plot dos resultados
figure(2);
subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada: impulso');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso: com ru�do');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');


%% 2) Usando entrada aleat�ria
b = 4 ; m =1;
u = prbs(N,b,m) - 0.5;
y = lsim(H, u, t);

% Estimando saida usando convolucao
Umat = zeros(N);
for i = 1:N
    Umat(i:end,i) = u(1:end-i+1)';
end
y_est = (Umat^-1)*y;

% Plot dos resultados
figure(3);
subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada: aleat�ria');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');

%% 2.1) Usando uma sa�da com ruido

y = y + 0.05*randn(length(y),1);
y_est = (Umat^-1)*y;


% Plot dos resultados
figure(4);
subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada: aleat�ria');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso: com ru�do');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');
