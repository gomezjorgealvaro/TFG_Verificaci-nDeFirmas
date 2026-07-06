% Metodo para tomar n individuos de la poblacion con muestreo proporcional
% al tamaño con reeemplazamiento
%X matriz de datos de la información auxiliar
%X matriz columna N*1
% n tamaño de la muestra
%M es el valor que garantiza que M>=max(X)
function muestra=lahiri(M,X,n)
N=size(X,1);
U=(1:1:N);
j=1;
muestra=zeros(n,1);
while (j<=n)
    %primer valor tomado aleatoriemente para escoger el individuo de la
    %población
ind=fix(unifrnd(1,N+1));
%segundo valor que se toma para decidir si el individuo es seleccionado en
%la muestra o no comparandolo con el tamaño del individuo ind

ind2=unifrnd(0,M);
if ind2<=X(ind,1)
muestra(j,1)=U(1,ind);
j=j+1;
end
end

