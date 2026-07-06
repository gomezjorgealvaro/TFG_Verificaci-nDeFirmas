function en = nEsperado(p, n1, c1, r1, n2)
% Tamano muestral esperado del plan doble (Seccion 2.3.4)
    q    = 1 - p;
    % Probabilidad de pasar a segunda etapa
    p_II = 0;
    for d1 = c1+1 : r1-1
        p_II = p_II + binopdf(d1, n1, q);
    end
    en = n1 + n2 * p_II;
end