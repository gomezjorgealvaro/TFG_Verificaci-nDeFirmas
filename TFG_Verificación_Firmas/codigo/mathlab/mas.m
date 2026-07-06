% Obtención de una m.a.s. de una población
% N: Tamaño de la población
% n: Tamaño de la muestra

function muestra=mas(N,n)
if n<=N
un=fix(unifrnd(1,N+1));
j=1;
muestra(1,1)=un;
while (j<n);
	un=fix(unifrnd(1,N+1));
if un~=muestra(1:length(muestra),1)
	j=j+1;
	muestra(j,1)=un;
end
end
end


