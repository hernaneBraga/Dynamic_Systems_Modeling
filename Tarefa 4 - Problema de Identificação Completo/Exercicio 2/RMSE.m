function rmse = RMSE(yreal, yhat)
    erro = yreal - yhat;
    rmse = sqrt(sum(erro.^2)) / sqrt(sum((yreal - mean(yreal)).^2));
end

