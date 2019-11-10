% Hernane Braga Pereira - UFMG Agosto 2019
% Exerc�cio 4.13 - Livro Introdu��o � Identifica��o de Sistemas 3.Ed

% Implementa��o das rotinas usadas s�o de autoria do autor do livro:
% (c) Luis Aguirre, 1999, 2011, 2012
% Dispon�veis em: https://www.researchgate.net/publication/303679484_Introducao_a_Identificacao_de_Sistemas

%%
clear
close all

N = 200; %num amostras
tb = 10; %intervalo entre bits
b = 6; %numero de bits

% Gera��o da PRBS
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
title('Autocorrela��o do sinal');
axis([-N N -0.5 1.5]);

