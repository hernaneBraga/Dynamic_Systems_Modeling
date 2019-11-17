function e = vetor_erro(y, S_R)
% y - Vetor de resposta
% S_R - relacao sinal/ruido em db
    N = length(y);
    sigma_y = std(y);
    sigma_e = sigma_y/(10^(S_R/20));
    e = normrnd(0, sigma_e, [1, N]);
    e = e';

end

