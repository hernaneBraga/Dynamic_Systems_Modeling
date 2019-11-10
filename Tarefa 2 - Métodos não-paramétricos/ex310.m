% Hernane Braga Pereira - UFMG Agosto 2019
% Exercício 3.10 - Livro Introdução à Identificação de Sistemas 3.Ed

% Implementação das rotinas usadas são de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Disponíveis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clear
close all

% Definicao da funcao de transferencia
H = tf(1, [1 0.2 0.8]);

% Parametros usados:
N = 15; % Número de amostras
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
title('Resposta ao impulso: com ruído');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');


%% 2) Usando entrada aleatória
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
title('Entrada: aleatória');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');

%% 2.1) Usando uma saída com ruido

y = y + 0.05*randn(length(y),1);
y_est = (Umat^-1)*y;


% Plot dos resultados
figure(4);
subplot(212)
plot(t, u);
xlabel('Tempo (s)')
ylabel('Amplitude');
title('Entrada: aleatória');

subplot(211)
plot(t, y, 'b');
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Resposta ao impulso: com ruído');
hold on
plot(t, y_est, 'ro');
legend('H(t)','He(t)');
