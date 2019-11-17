function y_sim_livre = validacao_livre(u, yreal, params)
%UNTITLED3 Summary of this function goes here
    N = length(yreal);
    
    y_sim_livre = zeros(N - 3, 1);
    y_sim_livre(1) = [yreal(2) yreal(1) 1 1] * params;
    y_sim_livre(2) = [y_sim_livre(1) yreal(2) 1 1] * params;
    
    for k = 3:(N-3)
        y_sim_livre(k) = [y_sim_livre(k-1) y_sim_livre(k-2) 1 1] * params;
    end


end

