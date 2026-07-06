function p = poisbinom_atleast(pis, kmin)
% POISBINOM_ATLEAST  Calcula P(K >= kmin), donde K es el numero de exitos
% de H ensayos de Bernoulli independientes con probabilidades distintas
% pis(1),...,pis(H) (distribucion de Poisson-Binomial), mediante
% convolucion iterativa de las distribuciones individuales.
%
% pis  : vector de probabilidades de exito (una por provincia, pi_h)
% kmin : numero minimo de exitos exigido (en este TFG, kmin = 5)

    H = length(pis);
    pmf = 1;  % P(K=0)=1 antes de incorporar ninguna provincia
    for h = 1:H
        pmf = conv(pmf, [1 - pis(h), pis(h)]);
    end
    % pmf(k+1) = P(K = k), para k = 0,...,H
    p = sum(pmf(kmin+1:end));
end