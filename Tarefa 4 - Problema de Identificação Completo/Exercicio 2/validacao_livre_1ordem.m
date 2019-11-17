function y_sim_livre = validacao_livre_1ordem(u, yreal, params)
%UNTITLED3 Summary of this function goes here
    N = length(yreal);
    
    y_sim_livre = zeros(N - 3, 1);
    y_sim_livre(1) = [yreal(1) u(1)] * params;
    y_sim_livre(2) = [y_sim_livre(1) u(2)] * params;
    
    for k = 3:(N-3)
        y_sim_livre(k) = [y_sim_livre(k-1) u(k)] * params;
    end


end