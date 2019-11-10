% Hernane Braga Pereira - UFMG Setembro 2019
% Exerc�cio 4.8 - Livro Introdu��o � Identifica��o de Sistemas 3.Ed

% Implementa��o das rotinas usadas s�o de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Dispon�veis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clc;
clear;
close all;

sd = 1;     % Desvio padr�o
med = 0;    % M�dia
N = 10000;  % N�mero de amostras

e = normrnd(med, sd, [1, N+3]); % Ru�do branco (gaussina de med e sd)
u = zeros(1, N-3);

% Calculando o valor do sinal u(k) pedido
for k=4:N+3
    u(k-3)= 0.9*e(k-1) + 0.8*e(k-2) + 0.7*e(k-3) + e(k);
end

% Plot do intervalo de confian�a
k = 5;
[t,ruu,l,B]=myccf2([u' u'],k,0,1,'k');
plot([0 5],[l l],'b-',[0 5],[-l -l],'b-');
axis([-0.5 5.5 -0.25 1.25]);
hold on

% Plot da autocorrela��o calculada
stem(t,ruu,'r');
hold off
xlabel('Atraso');
ylabel('r_u_u(k)');
title('Autocorrela��o de u(k)');


disp('Autocorrela��o m�xima B:');
disp(B)

disp('r_uu(k): ')
disp(ruu*B)

