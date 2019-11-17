function y_passo_frente = validacao_passo_frente(u2, yreal, ordem, params)
    N = length(yreal);
    y_passo_frente = [];
    for i = 1:N
        psi_p = []; % Vetor de regressores
        for j = 1:ordem
            if j >= i
                psi_p = [psi_p, yreal(i), u2(i)];
            else
                psi_p = [psi_p, yreal(i - j), u2(i - j)]; 
            end
        end
        y_passo_frente(i) = psi_p*params;
    end
end

