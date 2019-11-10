% Hernane Braga Pereira - UFMG Agosto 2019
% Exercício 4.18 - Livro Introdução à Identificação de Sistemas 3.Ed

% Implementação das rotinas usadas são de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Disponíveis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%

close all;
clear;
clc;

% Definição da função de transferência
num = 1;
den = [1000 1];
H = tf(num,den);
Ts = 1;

% Vetor com os valores de Tb pedidos
Tb = [1 100 1000 10000];
b = 12;
N = 10000;

% Plot da entrada e saída do sistema para cada Tb
for i=1:length(Tb)
    u = prbs(N, b, Tb(i));
    figure(i)
    lsim(H, u, 0:Ts:N-1, 'r-');
    xlim([0 N])
    title(['Resposta do sistema. Tb = ' num2str(Tb(i))]);
    xlabel('Tempo (s)')
end




