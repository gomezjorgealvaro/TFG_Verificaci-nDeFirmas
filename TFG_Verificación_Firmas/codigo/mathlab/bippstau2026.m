function [estau, estvtau]=bippstau2026(pk, nI, estauyi)
estau=(1/nI)*sum(estauyi./pk);
aux1 = (estauyi./pk).^2;
aux1 = sum(aux1);
aux2 = nI*estau*estau;
estvtau = (aux1-aux2)/(nI*(nI-1));
end

