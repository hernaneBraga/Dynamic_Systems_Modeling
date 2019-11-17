function y_sim_livre = validacao_livre_3ordem(u, yreal, params)
%UNTITLED3 Summary of this function goes here
    N = length(yreal);

    y_sim_livre = zeros(N - 3, 1);
    y_sim_livre(3) = [yreal(3) yreal(2) yreal(1) u(3) u(2) u(1)] * params;
    y_sim_livre(4) = [y_sim_livre(3) yreal(3) yreal(2) u(4) u(3) u(2)] * params;
    y_sim_livre(5) = [y_sim_livre(4) y_sim_livre(3) yreal(3) u(5) u(4) u(3)] * params;
    y_sim_livre(6) = [y_sim_livre(5) y_sim_livre(4) y_sim_livre(3) u(6) u(5) u(4)] * params;
    
    for k = 7:(N-3)
        y_sim_livre(k) = [y_sim_livre(k-1) y_sim_livre(k-2) y_sim_livre(k-3) u(k) u(k-1) u(k-2)] * params;
    end



end

