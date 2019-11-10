% Hernane Braga Pereira - UFMG Setembro 2019
% Exercício 4.8 - Livro Introdução à Identificação de Sistemas 3.Ed

% Implementação das rotinas usadas são de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Disponíveis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clc;
clear;
close all;

sd = 1;     % Desvio padrão
med = 0;    % Média
N = 10000;  % Número de amostras

e = normrnd(med, sd, [1, N+3]); % Ruído branco (gaussina de med e sd)
u = zeros(1, N-3);

% Calculando o valor do sinal u(k) pedido
for k=4:N+3
    u(k-3)= 0.9*e(k-1) + 0.8*e(k-2) + 0.7*e(k-3) + e(k);
end

% Plot do intervalo de confiança
k = 5;
[t,ruu,l,B]=myccf2([u' u'],k,0,1,'k');
plot([0 5],[l l],'b-',[0 5],[-l -l],'b-');
axis([-0.5 5.5 -0.25 1.25]);
hold on

% Plot da autocorrelação calculada
stem(t,ruu,'r');
hold off
xlabel('Atraso');
ylabel('r_u_u(k)');
title('Autocorrelação de u(k)');


disp('Autocorrelação máxima B:');
disp(B)

disp('r_uu(k): ')
disp(ruu*B)

