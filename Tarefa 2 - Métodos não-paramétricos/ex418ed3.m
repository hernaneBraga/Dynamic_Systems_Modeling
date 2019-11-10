% Hernane Braga Pereira - UFMG Agosto 2019
% Exerc�cio 4.18 - Livro Introdu��o � Identifica��o de Sistemas 3.Ed

% Implementa��o das rotinas usadas s�o de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Dispon�veis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%

close all;
clear;
clc;

% Defini��o da fun��o de transfer�ncia
num = 1;
den = [1000 1];
H = tf(num,den);
Ts = 1;

% Vetor com os valores de Tb pedidos
Tb = [1 100 1000 10000];
b = 12;
N = 10000;

% Plot da entrada e sa�da do sistema para cada Tb
for i=1:length(Tb)
    u = prbs(N, b, Tb(i));
    figure(i)
    lsim(H, u, 0:Ts:N-1, 'r-');
    xlim([0 N])
    title(['Resposta do sistema. Tb = ' num2str(Tb(i))]);
    xlabel('Tempo (s)')
end




