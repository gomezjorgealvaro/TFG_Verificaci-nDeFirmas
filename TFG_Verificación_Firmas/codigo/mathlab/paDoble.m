function pa = paDoble(p, n1, c1, r1, n2, c2)
% Probabilidad de aceptacion del plan doble (Proposicion 2.3.3)
    q  = 1 - p;
    % Aceptacion en primera etapa
    pa1 = binocdf(c1, n1, q);
    % Aceptacion en segunda etapa (condicionada)
    pa2 = 0;
    for d1 = c1+1 : r1-1
        p_d1   = binopdf(d1, n1, q);
        % D1+D2 <= c2 con D1=d1 => D2 <= c2-d1
        resto  = c2 - d1;
        if resto >= 0
            p_d2 = binocdf(resto, n2, q);
        else
            p_d2 = 0;
        end
        pa2 = pa2 + p_d1 * p_d2;
    end
    pa = pa1 + pa2;
end