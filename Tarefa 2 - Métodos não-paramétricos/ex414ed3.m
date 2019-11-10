% Hernane Braga Pereira - UFMG Agosto 2019
% Exercício 4.13 - Livro Introdução à Identificação de Sistemas 3.Ed

% Implementação das rotinas usadas são de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Disponíveis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clear
close all

N = 200; %num amostras
tb = 10; %intervalo entre bits
b = 6; %numero de bits

% Geração da PRBS
u=prbs(N*tb, b, tb);

figure(1)

subplot(211)
plot(u)
title('PRBS Tb = 10');
xlabel('Amostras')
axis([0 N -0.5 1.5]);

subplot(212)
[t,ruu,l,B]=myccf2(u',N*2,1,1,'k');
xlabel('Atraso')
title('Autocorrelação do sinal');
axis([-N N -0.5 1.5]);

