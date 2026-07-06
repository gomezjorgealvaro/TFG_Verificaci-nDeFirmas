function [estPest,estvPest, icPest] = estEstratificado(prov, X, Nh, nh, Wh, k)
idx_muestra = [];
H = length(Nh);
p_h = zeros(H,1);
var_h = zeros(H,1);
for h = 1:H
    idx_estrato = find(prov == h);
    N_h = length(idx_estrato);
    n_h = nh(h);
    if n_h > 0
        idx_sel = randsample(idx_estrato, n_h);
        idx_muestra = [idx_muestra; idx_sel];
        x_h = X(idx_sel);
        p_h(h) = mean(x_h);
        f_h = n_h / N_h;
        var_h(h) = (1 - f_h) * p_h(h) * (1 - p_h(h)) / (n_h - 1);
    else
        p_h(h) = NaN;
        var_h(h) = 0;
    end
end

Wh = Wh(:);        
p_h = p_h(:);        
var_h = var_h(:);

estPest = sum(Wh.*p_h);
estvPest = sum((Wh).^2.*var_h);
BPest = k*sqrt(estvPest);
icPest = [estPest-BPest estPest+BPest];
end