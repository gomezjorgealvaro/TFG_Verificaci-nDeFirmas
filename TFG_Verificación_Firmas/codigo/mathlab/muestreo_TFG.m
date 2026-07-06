
% Parámetros
N = 600000;          % Total de firmas presentadas
p_real = 0.85;        % Proporción real de firmas válidas (desconocida)
p0 = 0.8333;         % Umbral legal mínimo (500000/600000)
alpha = 0.05;        % Nivel de significación
B = 0.02;            % Error de estimación 
k = norminv(1 - alpha/2);   %  Valor crítico Aprox. 1.96

% Generación de la población (X=1 válida, X=0 inválida)
X = binornd(1, p_real, N, 1);   % Vector columna de N muestras Bernoulli
n_validas = sum(X);
Xp = n_validas / N;

% Resultados de la población
fprintf(' POBLACION SIMULADA DE FIRMAS ILP  \n');
fprintf(' Total firmas (N)            : %d\n', N);
fprintf(' Firmas válidas              : %d (%.4f)\n', n_validas, Xp);
fprintf(' Firmas inválidas            : %d (%.4f)\n', N-n_validas, 1-Xp);
fprintf(' Umbral legal (p0)           : %.4f\n', p0);
fprintf(' ¿Supera umbral?             : %s\n', string( Xp >= p0 ));


%MAS
% Cálculo del tamaño muestral para MAS (Muestreo Aleatorio Simple)
% Con p conservador = 0.5
n0 = (k*k * 0.25 * N) / (B^2 * (N-1));
n = ceil( n0 / (1 + n0/N) );

% Con p estimado de muestra previa = 0.85
n0_p = (k*k * 0.85 * 0.15 * N) / (B^2 * (N-1));
n_MASpre = ceil( n0_p / (1 + n0_p/N) );

fprintf('\n Tamaño muestral con p conservador (0.5): %d\n', n);
fprintf(' Tamaño muestral con p = 0.85           : %d\n', n_MASpre);

Rb = 1000;

contm = 0;
contmp = 0;

estP_acum  = 0;  icP_acum  = [0 0];
estPp_acum = 0;  icPp_acum = [0 0];

for b = 1:Rb
m1 = mas(N, n);
xm1 = X(m1, 1);

estP = mean(xm1);
f = n/N;
estvP = (1-f)*estP*(1-estP)/(n-1);
BP = k*sqrt(estvP);
icP = [estP-BP estP+BP];

mp1 = mas(N, n_MASpre);
xmp1 = X(mp1, 1);

estPp = mean(xmp1);
fp = n_MASpre/N;
estvPp = (1-fp)*estPp*(1-estPp)/(n_MASpre-1);
BPp = k*sqrt(estvPp);
icPp = [estPp-BPp estPp+BPp];

estP_acum  = estP_acum  + estP;
icP_acum   = icP_acum   + icP;
estPp_acum = estPp_acum + estPp;
icPp_acum  = icPp_acum  + icPp;

if p_real> estP-BP & p_real<estP+BP
    contm=contm+1;
end

if p_real> estPp-BPp & p_real<estPp+BPp
    contmp=contmp+1;
end

end

estP  = estP_acum  / Rb;
icP   = icP_acum   / Rb;
estPp = estPp_acum / Rb;
icPp  = icPp_acum  / Rb;
estvP  = (1-f)*estP*(1-estP)/(n-1);
estvPp = (1-fp)*estPp*(1-estPp)/(n_MASpre-1);

fprintf('\n Para la muestra con solución conservadora  %.4f\n', n);
fprintf('Estimador de P (promedio Monte Carlo, %d repl.):  %.4f\n', Rb, estP);
fprintf('Estimador de Varianza: %.4f\n', estvP);
fprintf('IC 95%% (promedio): [%.4f, %.4f]\n', icP);
se_H0 = sqrt( p0*(1-p0)/n * ( (N - n) / (N - 1) ) );
Z = (estP - p0) / se_H0;   % Estadístico de prueba (normal)
p_valor = normcdf(Z);      % Unilateral izquierda

fprintf('Estadístico Z = %.4f\n', Z);
fprintf('p-valor (unilateral) = %.6f\n', p_valor);
if p_valor < alpha
    fprintf('Rechazamos H0: la proporción de firmas válidas es significativamente inferior a %.4f\n', p0);
else
    fprintf('No se rechaza H0: no hay evidencia suficiente para afirmar que la proporción sea inferior a %.4f\n', p0);
end
fprintf('Con una confianza del  %.3f (%d de %d)\n', contm/Rb, contm, Rb);


fprintf('\n Para la muestra con proporción previa  %.4f\n', n_MASpre);
fprintf('Estimador de P (promedio Monte Carlo, %d repl.):  %.4f\n', Rb, estPp);
fprintf('IC 95%% (promedio): [%.4f, %.4f]\n', icPp);
se_H0p = sqrt( p0*(1-p0)/n_MASpre * ( (N - n_MASpre) / (N - 1) ) );
Zp = (estPp - p0) / se_H0p;
p_valorp = normcdf(Zp);
fprintf('Estadístico Z = %.4f\n', Zp);
fprintf('p-valor (unilateral) = %.6f\n', p_valorp);
if p_valorp < alpha
    fprintf('Rechazamos H0: la proporción de firmas válidas es significativamente inferior a %.4f\n', p0);
else
    fprintf('No se rechaza H0: no hay evidencia suficiente para afirmar que la proporción sea inferior a %.4f\n', p0);
end
fprintf('Con una confianza del  %.3f (%d de %d)\n', contmp/Rb, contmp, Rb);

%ESTRATIFICADO (9 PROVINCIAS)
% Datos recogidos de la Junta de Castilla y Leon a 1/1/2025 (más
% actualizados)
provincias = {"Ávila", "Burgos", "León", "Palencia", "Salamanca", ...
              "Segovia", "Soria", "Valladolid", "Zamora"};
datos = [159887 361556 446771 158687 328779 158470 90234 528841 165275];
Nh = round(datos/sum(datos)*N);
Wh = Nh/N;

dif = N-sum(Nh);
if dif > 0
    [~, idx] = max(Nh);
    Nh(idx) = Nh(idx) + dif;
elseif dif < 0
    [~, idx] = min(Nh);
    Nh(idx) = Nh(idx) + dif;
end

% Generar vector de provincias (repetir cada provincia según Nh)
Rb = 1000;
contest1 = 0;
contest2 = 0;

estPest1_acum = 0; estvPest1_acum = 0; icPest1_acum = [0 0];
estPest2_acum = 0; estvPest2_acum = 0; icPest2_acum = [0 0];

for b = 1:Rb
prov = [];
for h = 1:9
    prov = [prov; h * ones(Nh(h), 1)];
end
% Generar variable de firma válida (X) independiente de la provincia
 
X = binornd(1, p_real, N, 1);
p_real_actual = mean(X);
%Para Asignación proporcional usando el tamaño muestral de MAS calculamos 
%el tamaño muestral para cada estrato
nh = round(n*Nh/N);

dif = n - sum(nh);
if dif > 0
   [~, idx] = max(nh);
   nh(idx) = nh(idx) + dif;
elseif dif < 0
   [~, idx] = min(nh);
   nh(idx) = nh(idx) + dif;
end
 
[estPest1,estvPest1, icPest1] = estEstratificado(prov, X, Nh, nh, Wh, k);
if p_real_actual > icPest1(1) & p_real_actual < icPest1(2)
    contest1=contest1+1;
end
estPest1_acum  = estPest1_acum  + estPest1;
estvPest1_acum = estvPest1_acum + estvPest1;
icPest1_acum   = icPest1_acum   + icPest1;

%Para asignación óptima de Neyman

Sh = sqrt(p_real*(1-p_real)*Nh/(Nh-1));

num = Nh*Sh;
nh = round(n*num/sum(num));

dif = n - sum(nh);
if dif > 0
   [~, idx] = max(nh);
   nh(idx) = nh(idx) + dif;
elseif dif < 0
   [~, idx] = min(nh);
   nh(idx) = nh(idx) + dif;
end

[estPest2,estvPest2, icPest2] = estEstratificado(prov, X, Nh, nh, Wh, k);
if p_real_actual > icPest2(1) & p_real_actual < icPest2(2)
    contest2=contest2+1;
end
estPest2_acum  = estPest2_acum  + estPest2;
estvPest2_acum = estvPest2_acum + estvPest2;
icPest2_acum   = icPest2_acum   + icPest2;
end

estPest1  = estPest1_acum  / Rb;
estvPest1 = estvPest1_acum / Rb;
icPest1   = icPest1_acum   / Rb;
estPest2  = estPest2_acum  / Rb;
estvPest2 = estvPest2_acum / Rb;
icPest2   = icPest2_acum   / Rb;

fprintf('\n Muestreo estratificado (asignación proporcional)\n');
fprintf('estP (promedio Monte Carlo, %d repl.) = %.4f\n', Rb, estPest1);
fprintf('IC 95%% (promedio) = [%.4f, %.4f]\n', icPest1(1), icPest1(2));
fprintf('Varianza media estimada = %.6f\n', estvPest1);
fprintf('Con una confianza del  %.3f (%d de %d)\n', contest1/Rb, contest1, Rb)


fprintf('\n Muestreo estratificado (asignación óptima de Neyman)\n');
fprintf('estP (promedio Monte Carlo, %d repl.) = %.4f\n', Rb, estPest2);
fprintf('IC 95%% (promedio) = [%.4f, %.4f]\n', icPest2(1), icPest2(2));
fprintf('Varianza media estimada = %.6f\n', estvPest2);
fprintf('Con una confianza del  %.3f (%d de %d)\n', contest2/Rb, contest2, Rb)

%En estos casos de estratificado no hacemos que exista heterogeneidad entre
%estratos por lo que nuestro muestro apenas reduce la varianza sobre MAS,
%para ver la ventaja real de este tipo de muestreo, haremos que la
%proporción de firmas no sea similar en todas las povincias

p_h_inicial = [0.74; 0.94; 0.80; 0.95; 0.69; 0.90; 0.86; 0.97; 0.76]; 
Nh = Nh(:);               
Wh = Nh / sum(Nh);            
p_h_inicial = p_h_inicial(:);        
media_ponderada = sum(Nh .* p_h_inicial) / sum(Nh);
p_h = p_h_inicial * (p_real / media_ponderada);
p_h = min(0.99, max(0.01, p_h));

Rb = 1000;         
n_estrat = n;         

cobertura_prop = 0;
cobertura_neyman = 0;
var_prop_acum = 0;
var_neyman_acum = 0;
estP_prop_acum = 0;
estP_neyman_acum = 0;

ic_prop_acum   = [0 0];
ic_neyman_acum = [0 0];

for b = 1:Rb
    prov_het = [];
    for h = 1:9
        prov_het = [prov_het; h * ones(Nh(h), 1)];
    end
    
    X_het = zeros(N, 1);
    for h = 1:9
        idx = (prov_het == h);
        n_h_estrato = sum(idx);
        X_het(idx) = binornd(1, p_h(h), n_h_estrato, 1);
    end
    p_real_actual_het = mean(X_het);  
    
    %  Asignación proporcional 
    nh_prop = round(n_estrat * Nh / N);
    diff_prop = n_estrat - sum(nh_prop);
    if diff_prop > 0
        [~, idx_max] = max(nh_prop);
        nh_prop(idx_max) = nh_prop(idx_max) + diff_prop;
    elseif diff_prop < 0
        [~, idx_min] = min(nh_prop);
        nh_prop(idx_min) = nh_prop(idx_min) + diff_prop;
    end
    
    [estP_prop, var_prop, ic_prop] = estEstratificado(prov_het, X_het, Nh, nh_prop, Wh, k);
    if p_real_actual_het >= ic_prop(1) && p_real_actual_het <= ic_prop(2)
        cobertura_prop = cobertura_prop + 1;
    end
    var_prop_acum = var_prop_acum + var_prop;
    estP_prop_acum = estP_prop_acum + estP_prop;
    ic_prop_acum = ic_prop_acum + ic_prop;
    
    % Asignación óptima de Neyman 
    Sh = sqrt(p_h .* (1 - p_h) .* Nh ./ (Nh - 1));  
    num = Nh .* Sh;
    nh_neyman = round(n_estrat * num / sum(num));
    diff_neyman = n_estrat - sum(nh_neyman);
    if diff_neyman > 0
        [~, idx_max] = max(nh_neyman);
        nh_neyman(idx_max) = nh_neyman(idx_max) + diff_neyman;
    elseif diff_neyman < 0
        [~, idx_min] = min(nh_neyman);
        nh_neyman(idx_min) = nh_neyman(idx_min) + diff_neyman;
    end
    
    [estP_neyman, var_neyman, ic_neyman] = estEstratificado(prov_het, X_het, Nh, nh_neyman, Wh, k);
    if p_real_actual_het >= ic_neyman(1) && p_real_actual_het <= ic_neyman(2)
        cobertura_neyman = cobertura_neyman + 1;
    end
    var_neyman_acum = var_neyman_acum + var_neyman;
    estP_neyman_acum = estP_neyman_acum + estP_neyman;
    ic_neyman_acum = ic_neyman_acum + ic_neyman;
end

ic_prop   = ic_prop_acum   / Rb;
ic_neyman = ic_neyman_acum / Rb;

% Mostrar resultados
fprintf('\n--- Resultados con heterogeneidad real ---\n');
fprintf('Tamaño muestral total utilizado: %d\n', n_estrat);

fprintf('\nESTRATIFICADO PROPORCIONAL:\n');
fprintf('  Con una confianza empírica: %.3f (%d/%d)\n', cobertura_prop/Rb, cobertura_prop, Rb);
fprintf('  Varianza media estimada: %.6f\n', var_prop_acum / Rb);
fprintf('  Error estándar medio (raíz(var)): %.6f\n', sqrt(var_prop_acum / Rb));
fprintf('  IC 95%% = [%.4f, %.4f]\n', ic_prop(1), ic_prop(2));
fprintf('  Estimador proporcion de p: %.4f\n', estP_prop_acum / Rb);

fprintf('\nESTRATIFICADO ÓPTIMO (NEΥMAN):\n');
fprintf('  Con una confianza empírica: %.3f (%d/%d)\n', cobertura_neyman/Rb, cobertura_neyman, Rb);
fprintf('  Varianza media estimada: %.6f\n', var_neyman_acum / Rb);
fprintf('  Error estándar medio (raíz(var)): %.6f\n', sqrt(var_neyman_acum / Rb));
fprintf('  IC 95%% = [%.4f, %.4f]\n', ic_neyman(1), ic_neyman(2));
fprintf('  Estimador proporcion de p: %.4f\n', estP_neyman_acum / Rb);

reduccion_var = (var_prop_acum - var_neyman_acum) / var_prop_acum * 100;
fprintf('\nLa asignación de Neyman reduce la varianza en un %.2f%% respecto a la proporcional.\n', reduccion_var);


%MUESTREO POR CONGLOMERADOS BIETÁPICO CON PPT
municipios_totales = 2248;   % Número real aproximado de municipios CyL

% Costes
c1 = 150;    % coste fijo por desplazamiento al municipio
c2 = 2;      % coste por inspección de firma


rho = 0.03;%Correlación entre municipios
%Simulamos tamaños muestrales
% Distribución muy asimétrica: muchos municipios pequeños
mu_ln = 4.5;
sigma_ln = 1.2;
x = round(lognrnd(mu_ln, sigma_ln, municipios_totales, 1));
% Evitar municipios demasiado pequeños
x(x < 20) = 20;
% Total de firmas simuladas
N_total = sum(x);

fprintf('\n----- POBLACIÓN MUNICIPAL SIMULADA -----\n');
fprintf('Municipios totales              : %d\n', municipios_totales);
fprintf('Total firmas simuladas          : %d\n', N_total);
fprintf('Media firmas/municipio          : %.2f\n', mean(x));
fprintf('Máximo tamaño municipal         : %d\n', max(x));
fprintf('Mínimo tamaño municipal         : %d\n', min(x));

%Calculamos tamaó óptimo según costes
k_opt = round(sqrt(c1*(1-rho)/(c2*rho)));
DEFF_cost = 1 + (k_opt - 1)*rho;

n_cong_necesario_cost = ceil(n * DEFF_cost);
nI_cost    = ceil(n_cong_necesario_cost / k_opt);
n_cong_cost = nI_cost * k_opt;

fprintf('\n----- DISEÑO ÓPTIMO -----\n');
fprintf('rho intraclase                  : %.3f\n', rho);
fprintf('k* óptimo                       : %d\n', k_opt);
fprintf('DEFF                            : %.3f\n', DEFF_cost);
fprintf('Municipios seleccionados        : %d\n', nI_cost);
fprintf('Tamaño muestral total           : %d\n', n_cong_cost);

%Generamos proporciones por municipio
% Variabilidad entre municipios
sigma_mun = sqrt(rho * p_real * (1-p_real));
p_municipios = p_real + sigma_mun*randn(municipios_totales,1);

% Acotar probabilidades
p_municipios = max(0.01, min(0.99, p_municipios));

%Realizamos una iteración de Monte-carlo para probar laa confianza
Rb_cong = 1000;
cont_cong_cost = 0;
estPscong_cost = zeros(Rb_cong,1);
Bscong_cost    = zeros(Rb_cong,1);
icPcong_cost_acum = [0 0];

for b = 1:Rb_cong
    %Primera etapa con muestreo proporcional al tamaño
    mun_sel = lahiri(max(x), x, nI_cost);
    % Tamaños de los conglomerados seleccionados
    xis = x(mun_sel);
    % Probabilidades de selección aproximadas
    pis = xis / N_total;
    %Segunda etapa MAS en cada municipio
    estauyi = zeros(nI_cost,1);
    for i = 1:nI_cost
        % Conglomerado seleccionado
        h = mun_sel(i);
        % Tamaño poblacional del conglomerado
        % CORREGIDO: se llamaba "Nh", igual que el vector provincial de
        % 9 elementos usado en los bloques de estratificado, lo que
        % sobrescribia esa variable global con un escalar. Se renombra
        % para evitar el shadowing.
        Nh_mun = x(h);
        % Tamaño muestral dentro del conglomerado
        nh = min(k_opt, Nh_mun);
        % Generar población del conglomerado
        firmas_mun = binornd(1, p_municipios(h), Nh_mun, 1);
        % MAS sin reemplazo
        sel = mas(Nh_mun, nh);
        muestra_h = firmas_mun(sel);
        % Estimador del total del conglomerado
        estauyi(i) = Nh_mun * mean(muestra_h);
    end
    [estTau, estvTau] = bippstau2026(pis, nI_cost, estauyi);
    % Estimador de proporción poblacional
    estPcong = estTau / N_total;
    % Varianza del estimador de proporción
    estvPcong = estvTau / (N_total^2);
    % IC 95%
    t_crit = tinv(1 - alpha/2, nI_cost - 1);
    Bcong = t_crit * sqrt(estvPcong);
    icPcong = [estPcong - Bcong, estPcong + Bcong];
    % Guardar resultados
    estPscong_cost(b) = estPcong;
    Bscong_cost(b) = Bcong;
    icPcong_cost_acum = icPcong_cost_acum + icPcong;
    % Cobertura
    if p_real >= icPcong(1) && p_real <= icPcong(2)
        cont_cong_cost = cont_cong_cost + 1;
    end
end

icPcong_cost = icPcong_cost_acum / Rb_cong;

fprintf('\n----- RESULTADOS MUESTREO POR CONGLOMERADOS CON COSTES -----\n');
fprintf('Estimación proporción de p           : %.4f\n', mean(estPscong_cost));
fprintf('  IC 95%% = [%.4f, %.4f]\n', icPcong_cost(1), icPcong_cost(2));
fprintf('Sesgo empírico                  : %.6f\n',mean(estPscong_cost) - p_real);
fprintf('Error estándar medio            : %.6f\n', mean(Bscong_cost)/t_crit);
fprintf('Amplitud media IC95%%           : %.6f\n',mean(2*Bscong_cost));
fprintf('Con una confianza empírica IC95%%        : %.3f (%d/%d)\n', ...
    cont_cong_cost/Rb_cong, cont_cong_cost, Rb_cong);

fprintf('DEFF teórico                    : %.3f\n', DEFF_cost);


%Solución sin costes

rho = 0.03;
DEFF_max = 1.3;  % tolerancia máxima de pérdida de eficiencia frente al MAS

k_max_admisible = 1 + (DEFF_max - 1)/rho;
k_fijo = floor(k_max_admisible);  % redondeando hacia abajo para no superar la cota

DEFF_err = 1 + (k_fijo - 1)*rho;
n_cong_necesario_err = ceil(n * DEFF_err);
nI_err = ceil(n_cong_necesario_err / k_fijo);


%Generamos proporciones por municipio
% Variabilidad entre municipios
sigma_mun = sqrt(rho * p_real * (1-p_real));
p_municipios = p_real + sigma_mun*randn(municipios_totales,1);

% Acotar probabilidades
p_municipios = max(0.01, min(0.99, p_municipios));

%Realizamos una iteración de Monte-carlo para probar laa confianza
Rb_cong = 1000;
cont_cong_err = 0;
estPscong_err = zeros(Rb_cong,1);
Bscong_err    = zeros(Rb_cong,1);
icPcong_err_acum = [0 0];

for b = 1:Rb_cong
    %Primera etapa con muestreo proporcional al tamaño
    mun_sel = lahiri(max(x), x, nI_err);
    % Tamaños de los conglomerados seleccionados
    xis = x(mun_sel);
    % Probabilidades de selección aproximadas
    pis = xis / N_total;
    %Segunda etapa MAS en cada municipio
    estauyi = zeros(nI_err,1);
    for i = 1:nI_err
        % Conglomerado seleccionado
        h = mun_sel(i);
        % Tamaño poblacional del conglomerado
        Nh_mun = x(h);
        % Tamaño muestral dentro del conglomerado
        nh = min(k_fijo, Nh_mun);
        % Generar población del conglomerado
        firmas_mun = binornd(1, p_municipios(h), Nh_mun, 1);
        % MAS sin reemplazo
        sel = mas(Nh_mun, nh);
        muestra_h = firmas_mun(sel);
        % Estimador del total del conglomerado
        estauyi(i) = Nh_mun * mean(muestra_h);
    end
    [estTau, estvTau] = bippstau2026(pis, nI_err, estauyi);
    % Estimador de proporción poblacional
    estPcong = estTau / N_total;
    % Varianza del estimador de proporción
    estvPcong = estvTau / (N_total^2);
    % IC 95%
    t_crit = tinv(1 - alpha/2, nI_err - 1);
    Bcong = t_crit * sqrt(estvPcong);
    icPcong = [estPcong - Bcong, estPcong + Bcong];
    % Guardar resultados
    estPscong_err(b) = estPcong;
    Bscong_err(b) = Bcong;
    icPcong_err_acum = icPcong_err_acum + icPcong;
    % Cobertura
    if p_real >= icPcong(1) && p_real <= icPcong(2)
        cont_cong_err = cont_cong_err + 1;
    end
end

icPcong_err = icPcong_err_acum / Rb_cong;

%En este caso fijamos una cota máxima de DEFFF, para no perder un
%porcentaje superior con respecto a MAS
fprintf('\n----- RESULTADOS MUESTREO POR CONGLOMERADOS SIN COSTES -----\n');
fprintf('Estimación proporción de p           : %.4f\n', mean(estPscong_err));
fprintf('  IC 95%% = [%.4f, %.4f]\n', icPcong_err(1), icPcong_err(2));
fprintf('Sesgo empírico                  : %.6f\n',mean(estPscong_err) - p_real);
fprintf('Error estándar medio            : %.6f\n', mean(Bscong_err)/t_crit);
fprintf('Amplitud media IC95%%           : %.6f\n',mean(2*Bscong_err));
fprintf('Con una confianza empírica IC95%%        : %.3f (%d/%d)\n', ...
    cont_cong_err/Rb_cong, cont_cong_err, Rb_cong);

fprintf('DEFF teórico                    : %.3f\n', DEFF_err);


%Para el muestreo secuencial de Wald haremos 4 simulaciones para comprobar
%todos los casos

%MUESTREO SECUENCIAL DE WALD (PRPS) 
%Parametros del contraste
p1_s    = p_real;   % H0: p >= p1 (ILP alcanza el umbral)
p2_s    = 0.80;     % H1: p  = p2 (ILP no alcanza el umbral)
alpha_s = 0.05;
beta_s  = 0.10;

%Umbrales de decision en escala logaritmica
A_log = log(beta_s  / (1 - alpha_s));    % Umbral inferior: aceptar H0
B_log = log((1 - beta_s) / alpha_s);     % Umbral superior: rechazar H0

%Incremento de log-verosimilitud por cada firma inspeccionada
lr1 = log(p2_s / p1_s);           % Contribucion si X_i = 1 (valida)
lr0 = log((1-p2_s) / (1-p1_s));   % Contribucion si X_i = 0 (invalida)

fprintf('\n ------MUESTREO SECUENCIAL DE WALD------- \n')
fprintf(' H0: p >= %.4f   vs   H1: p = %.4f\n', p1_s, p2_s);
fprintf(' alpha = %.2f,  beta = %.2f\n', alpha_s, beta_s);
fprintf(' Umbral inferior A = %.4f  (log: %.4f)\n', exp(A_log), A_log);
fprintf(' Umbral superior B = %.4f  (log: %.4f)\n', exp(B_log), B_log);

[dec_unica, n_par, tray] = prps(X, A_log, B_log, lr1, lr0);

if dec_unica == 1
    str_dec = 'Aceptar H0  ->  ILP alcanza el umbral de validez';
elseif dec_unica == -1
    str_dec = 'Rechazar H0  ->  ILP no alcanza el umbral de validez';
else
    str_dec = 'Indeciso (maximo de iteraciones alcanzado)';
end

fprintf('\n Ejecucion sobre poblacion (p_real = %.4f):\n', p_real);
fprintf(' Decision        : %s\n', str_dec);
fprintf(' Parada en n     : %d\n', n_par);
fprintf(' Ahorro vs MAS   : %d firmas (%.1f%%)\n', ...
        n - n_par, (n - n_par)/n * 100);

% Simulacion 1: distribucion de n_parada bajo p_real 
Rb_seq = 1000;
n_paradas  = zeros(Rb_seq, 1);
decs_seq   = zeros(Rb_seq, 1);

for b = 1:Rb_seq
    X_b = binornd(1, p_real, N, 1);
    [decs_seq(b), n_paradas(b), ~] = prps(X_b, A_log, B_log, lr1, lr0);
end

E_n_seq    = mean(n_paradas);
med_n_seq  = median(n_paradas);
pct_acepta = sum(decs_seq ==  1) / Rb_seq;
pct_rechaza= sum(decs_seq == -1) / Rb_seq;

fprintf('\n Simulacion bajo p_real = %.4f  (%d repl.):\n', p_real, Rb_seq);
fprintf(' E[n] empirico            : %.1f\n', E_n_seq);
fprintf(' Mediana de n             : %.1f\n', med_n_seq);
fprintf(' Ahorro medio vs MAS      : %.1f firmas (%.1f%%)\n', ...
        n_MASpre - E_n_seq, (n_MASpre - E_n_seq)/n_MASpre * 100);
fprintf(' %% Acepta H0              : %.1f%%\n', pct_acepta  * 100);
fprintf(' %% Rechaza H0             : %.1f%%\n', pct_rechaza * 100);

%Simulacion 2: verificacion del error tipo I rechazar H0 cuando es H0(alpha)
Rb_H0   = 1000;
rec_H0  = 0;
n_H0    = zeros(Rb_H0, 1);

for b = 1:Rb_H0
    X_b = binornd(1, p1_s, N, 1);
    [d, n_H0(b), ~] = prps(X_b, A_log, B_log, lr1, lr0);
    if d == -1
        rec_H0 = rec_H0 + 1;
    end
end

fprintf('\n Verificacion error tipo I  (p = p1 = %.4f, %d repl.):\n', p1_s, Rb_H0);
fprintf(' P(rechazar H0 | H0 cierta) = %.3f  (nominal alpha = %.2f)\n', ...
        rec_H0/Rb_H0, alpha_s);
fprintf(' E[n] bajo H0               = %.1f\n', mean(n_H0));

%Simulacion 3: verificacion del error tipo II aceptar H0 cuando es H1(beta)
Rb_H1   = 1000;
acep_H1 = 0;
n_H1    = zeros(Rb_H1, 1);

for b = 1:Rb_H1
    X_b = binornd(1, p2_s, N, 1);
    [d, n_H1(b), ~] = prps(X_b, A_log, B_log, lr1, lr0);
    if d == 1
        acep_H1 = acep_H1 + 1;
    end
end

fprintf('\n Verificacion error tipo II (p = p2 = %.4f, %d repl.):\n', p2_s, Rb_H1);
fprintf(' P(aceptar H0 | H1 cierta)  = %.3f  (nominal beta = %.2f)\n', ...
        acep_H1/Rb_H1, beta_s);
fprintf(' E[n] bajo H1               = %.1f\n', mean(n_H1));

%Simulacion 4: E[n] en funcion de p (curva de eficiencia) cuando p se aleja
%de p1
p_vals  = linspace(0.75, 1.00, 20);
E_n_p   = zeros(size(p_vals));
Rb_curv = 1000;

for i = 1:length(p_vals)
    n_tmp = zeros(Rb_curv, 1);
    for b = 1:Rb_curv
        X_tmp = binornd(1, p_vals(i), N, 1);
        [~, n_tmp(b), ~] = prps(X_tmp, A_log, B_log, lr1, lr0);
    end
    E_n_p(i) = mean(n_tmp);
end

% Grafico: E[n] en funcion de p 

plot(p_vals, E_n_p, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, ...
     'MarkerFaceColor','r');
hold on;
yline(n, 'b--', 'LineWidth', 2.0, ...
      'Label', sprintf('n MAS = %d', n), ...
      'LabelHorizontalAlignment','right');
xline(p1_s, 'k:', 'LineWidth', 1.5, ...
      'Label', sprintf('p_1 = %.2f', p1_s));
xline(p2_s, 'k:', 'LineWidth', 1.5, ...
      'Label', sprintf('p_2 = %.2f', p2_s));
xlabel('p real');
ylabel('E[n]');
title(sprintf('E[n] del PRPS en funcion de p (%d sim. por punto)', Rb_curv));
legend('E[n] PRPS', sprintf('n MAS = %d', n), ...
       'Location','north', 'FontSize', 8);
grid on;
ylim([0, n*1.3]);


%PLANES DE MUESTREO DE ACEPTACION
% Parametros del plan
p1_plan  = 0.85;   % NCA: probabilidad aceptacion >= 1-alpha_pl
p2_plan  = 0.80;   % NCL: probabilidad aceptacion <= beta_pl
alpha_pl = 0.05;
beta_pl  = 0.10;

fprintf('\n\n\n');
fprintf('PLANES DE MUESTREO DE ACEPTACION\n');
fprintf(' NCA p1=%.2f,  NCL p2=%.2f\n', p1_plan, p2_plan);
fprintf(' alpha=%.2f,   beta=%.2f\n', alpha_pl, beta_pl);
 
%Funciones de probabilidad de aceptacion
% Plan simple bajo aproximacion binomial
pa_simple = @(n_pl, c_pl, p) binocdf(c_pl, n_pl, 1-p);
 
%Plan simple bajo distribucion hipergeometrica exacta
pa_hiper = @(N_pop, n_pl, c_pl, p) hygecdf(c_pl, N_pop, ...
    round(N_pop*(1-p)), n_pl);
 
%Busqueda del plan simple optimo
plan_encontrado = false;
n_opt = NaN; c_opt = NaN;
 
for n_pl = 10:4000
    for c_pl = 0:n_pl-1
        pa1 = pa_simple(n_pl, c_pl, p1_plan);
        pa2 = pa_simple(n_pl, c_pl, p2_plan);
        if pa1 >= 1-alpha_pl && pa2 <= beta_pl
            n_opt = n_pl;
            c_opt = c_pl;
            plan_encontrado = true;
            break;
        end
    end
    if plan_encontrado
        break;
    end
end
 
fprintf('\n\n')
fprintf('-----PLAN SIMPLE DE ACEPTACION OPTIMO-------\n');
fprintf(' Plan optimo: (n* = %d,  c* = %d)\n', n_opt, c_opt);
fprintf(' Pa(p1=%.2f) = %.4f  (requerido >= %.2f)\n', ...
    p1_plan, pa_simple(n_opt, c_opt, p1_plan), 1-alpha_pl);
fprintf(' Pa(p2=%.2f) = %.4f  (requerido <= %.2f)\n', ...
    p2_plan, pa_simple(n_opt, c_opt, p2_plan), beta_pl);
 
%Plan doble
n1_d = floor(n_opt / 2);
c1_d = max(0, floor(c_opt * 0.6));
r1_d = c1_d + 3;
n2_d = n_opt - n1_d;
c2_d = c_opt;
r2_d = c2_d + 1;
 
% Probabilidad de aceptacion del plan doble
pa_doble = @(p) paDoble(p, n1_d, c1_d, r1_d, n2_d, c2_d);
 
% Tamano muestral esperado del plan doble
n_esperado = @(p) nEsperado(p, n1_d, c1_d, r1_d, n2_d);
 
fprintf('\n Plan doble: (n1=%d, c1=%d, r1=%d, n2=%d, c2=%d, r2=%d)\n', n1_d, c1_d, r1_d, n2_d, c2_d, r2_d);
fprintf(' Pa(p1) doble = %.4f\n', pa_doble(p1_plan));
fprintf(' Pa(p2) doble = %.4f\n', pa_doble(p2_plan));
fprintf(' E[n] doble en p1 = %.1f\n', n_esperado(p1_plan));
 
%Curvas OC
p_v   = linspace(0.70, 1.00, 300);
oc_s  = arrayfun(@(p) pa_simple(n_opt, c_opt, p),  p_v);
oc_d  = arrayfun(@(p) pa_doble(p),                  p_v);
oc_h  = arrayfun(@(p) pa_hiper(N, n_opt, c_opt, p), p_v);
ne_v  = arrayfun(@(p) n_esperado(p),                 p_v);
 
%Graficos 
figure('Name','Planes de Aceptacion','NumberTitle','off','Position',[100 100 1100 480]);
 
%Panel izquierdo: Curvas OC
subplot(1, 2, 1);
plot(p_v, oc_s, 'b-',  'LineWidth', 2.5, ...
     'DisplayName', sprintf('Simple (n=%d, c=%d)', n_opt, c_opt));
hold on;
plot(p_v, oc_d, 'r--', 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Doble (n1=%d,c1=%d,r1=%d)', n1_d, c1_d, r1_d));
plot(p_v, oc_h, 'g:',  'LineWidth', 2.0, ...
     'DisplayName', 'Simple - Hipergeometrica');
xline(p1_plan, 'Color', [0, 0, 0.5],   'LineStyle',':', 'LineWidth',1.5, ...
      'Label', sprintf('NCA p1=%.2f', p1_plan), ...
      'LabelVerticalAlignment','bottom');
xline(p2_plan, 'Color', [1, 0.5, 0],   'LineStyle',':', 'LineWidth',1.5, ...
      'Label', sprintf('NCL p2=%.2f', p2_plan), ...
      'LabelVerticalAlignment','bottom');
yline(1-alpha_pl, 'Color',[0.5 0.5 0.5], 'LineStyle','--', 'LineWidth',1, ...
      'Label', sprintf('1-alpha=%.2f', 1-alpha_pl), ...
      'LabelHorizontalAlignment','right');
yline(beta_pl, 'Color',[1 0.6 0.6], 'LineStyle','--', 'LineWidth',1, ...
      'Label', sprintf('beta=%.2f', beta_pl), ...
      'LabelHorizontalAlignment','right');
xlabel('Proporcion de firmas validas (p)');
ylabel('Pa(p)');
title('Curvas Caracteristicas de Operacion (OC)');
legend('Location','southeast', 'FontSize', 8);
grid on;
xlim([0.70 1.00]);
ylim([0 1.05]);
 
% Panel derecho: E[n] plan doble vs simple
subplot(1, 2, 2);
plot(p_v, ne_v, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Plan doble E[n]');
hold on;
yline(n_opt, 'b--', 'LineWidth', 2.0, ...
      'Label', sprintf('Plan simple n=%d', n_opt), ...
      'LabelHorizontalAlignment','right');
yline(n1_d, 'g:', 'LineWidth', 1.5, ...
      'Label', sprintf('n1=%d (1a etapa)', n1_d), ...
      'LabelHorizontalAlignment','right');
xline(p1_plan, 'Color', [0, 0, 0.5], 'LineStyle', ':', 'LineWidth', 1.5);
xline(p2_plan, 'Color', [1, 0.5, 0], 'LineStyle', ':', 'LineWidth', 1.5);
xlabel('p');
ylabel('E[n]');
title('E[n] - Plan doble vs simple');
legend('Location','north', 'FontSize', 9);
grid on;
xlim([0.70 1.00]);
 
sgtitle('Planes de Muestreo de Aceptacion — ILP Castilla y Leon', ...
        'FontSize', 12, 'FontWeight', 'bold');
 
exportgraphics(gcf, 'mi_figura.png', 'Resolution', 300);


coste_MAS   = 500 + 2 * n;
coste_est   = 500 + 150*9 + 2 * n;
coste_cong  = 500 + 150*120 + 2 * 5880;
coste_seq   = round(500 + 2 * E_n_seq);
coste_acep  = 500 + 2 * n_opt;                      
n_esperado_p1 = n_esperado(p1_plan);                
coste_acep_doble = 500 + 2 * n_esperado_p1;

fprintf('\n\n\n');
fprintf('         TABLA COMPARATIVA DE DISENOS\n');
fprintf('   N=%d, p=%.2f, epsilon=%.2f, alpha=%.2f\n', N, p_real, B, alpha);
fprintf('\n');
fprintf('%-22s %8s %10s %12s\n', 'Diseno', 'n', 'Coste(eur)', 'Complejidad');
fprintf('%s\n', repmat('-', 1, 57));
fprintf('%-22s %8d %10d %12s\n', 'MAS',             n,       coste_MAS,  'Baja');
fprintf('%-22s %8d %10d %12s\n', 'Estratificado',   n,       coste_est,  'Media');
fprintf('%-22s %8d %10d %12s\n', 'Conglomerados',   n_cong,  coste_cong, 'Media');
fprintf('%-22s %8.0f %10d %12s\n','Secuencial PRPS', E_n_seq, coste_seq,  'Alta');
fprintf('%-22s %8d %10d %12s\n','Aceptacion simple', n_opt,   coste_acep, 'Baja');
fprintf('%-22s %8.0f %10d %12s\n','Aceptacion doble', n_esperado_p1, coste_acep_doble, 'Media');
fprintf('%s\n', repmat('-', 1, 57));


var_MAS = estvPp;
var_PROP = var_prop_acum / Rb;
var_NEYMAN = var_neyman_acum / Rb;
var_CONG = mean(Bscong.^2)/(k^2);

ER_PROP = var_MAS / var_PROP;
ER_NEYMAN = var_MAS / var_NEYMAN;
ER_CONG = var_MAS / var_CONG;

fprintf('\n\n\n');
fprintf('Eficiencia Relativa respecto a MAS\n');
fprintf('Estratificado proporcional : %.3f\n', ER_PROP);
fprintf('Estratificado Neyman       : %.3f\n', ER_NEYMAN);
fprintf('Conglomerados              : %.3f\n', ER_CONG);

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('   SIMULACION TERRITORIAL - CASTILLA Y LEON (Apartado 4.3)\n');
fprintf('================================================================\n');

%Censo real y umbral legal
N_censo  = 2090634;          % Censo electoral autonomico real
p0_terr  = 0.0075;           % Umbral legal: 0,75% (global y por provincia)

% Reutilizamos la distribucion poblacional por provincia (variable 'datos'
% ya definida en el bloque de muestreo estratificado)
Nh_censo = round(datos/sum(datos)*N_censo);
dif = N_censo - sum(Nh_censo);
if dif > 0
    [~, idx] = max(Nh_censo);
    Nh_censo(idx) = Nh_censo(idx) + dif;
elseif dif < 0
    [~, idx] = min(Nh_censo);
    Nh_censo(idx) = Nh_censo(idx) + dif;
end
Wh_censo = Nh_censo / N_censo;

fprintf('\nCenso electoral autonomico (N)        : %d\n', N_censo);
fprintf('Umbral legal global y provincial (p0) : %.4f (%.2f%%)\n', p0_terr, p0_terr*100);

% Asignacion muestral (proporcional, presupuesto fijo)
n_total_terr = 90000;   % Presupuesto muestral total para la auditoria
nh_terr = round(n_total_terr * Nh_censo / N_censo);
dif = n_total_terr - sum(nh_terr);
if dif > 0
    [~, idx] = max(nh_terr);
    nh_terr(idx) = nh_terr(idx) + dif;
elseif dif < 0
    [~, idx] = min(nh_terr);
    nh_terr(idx) = nh_terr(idx) + dif;
end

fprintf('Tamano muestral total (asig. proporcional) : %d\n', sum(nh_terr));
fprintf('Tamano muestral bajo asignacion uniforme   : %d (n/9)\n', round(n_total_terr/9));

% Escenarios de proporciones provinciales (ph)
% Orden: Avila, Burgos, Leon, Palencia, Salamanca, Segovia, Soria, Valladolid, Zamora
ph_homogeneo     = [0.0120 0.0125 0.0130 0.0118 0.0122 0.0115 0.0128 0.0135 0.0119];
ph_concentracion = [0.0040 0.0200 0.0210 0.0035 0.0195 0.0038 0.0030 0.0220 0.0042];
ph_frontera      = [0.0078 0.0072 0.0080 0.0070 0.0076 0.0074 0.0079 0.0073 0.0077];

escenarios = {'Homogeneo', 'Concentracion', 'Frontera'};
ph_mat     = [ph_homogeneo; ph_concentracion; ph_frontera];

% Simulacion Monte Carlo por escenario
Rb_terr = 2000;
resultados_terr = struct();
K_vals_all = zeros(Rb_terr, length(escenarios));   % para el grafico final

for e = 1:length(escenarios)

    ph_real = ph_mat(e, :);
    p_real_terr = sum(Wh_censo .* ph_real);   % proporcion global real ponderada

    aceptaciones = zeros(Rb_terr,1);
    K_vals       = zeros(Rb_terr,1);
    supera_h     = zeros(Rb_terr,9);   % 1 si la provincia h supera p0 en esa replica

    for b = 1:Rb_terr

        dh = binornd(nh_terr, ph_real);        % nº de firmas validas detectadas/provincia
        phat_h = dh ./ nh_terr;

        % Decision provincial (Definicion 2.5.1, condicion ii)
        decision_h = phat_h >= p0_terr;
        supera_h(b,:) = decision_h;
        K = sum(decision_h);
        K_vals(b) = K;

        % Estimador global estratificado (Ecuacion 2.13)
        phat_global = sum(Wh_censo .* phat_h);

        % Regla de aceptacion territorial (Definicion 2.5.1)
        aceptaciones(b) = (phat_global >= p0_terr) && (K >= 5);
    end

    K_vals_all(:,e) = K_vals;

    % Frecuencias empiricas pi_h: con que frecuencia cada provincia supera su umbral
    pi_h_emp = mean(supera_h, 1);

    % P(K>=5) teorica via Poisson-Binomial usando las pi_h empiricas
    pK5_teorica = poisbinom_atleast(pi_h_emp, 5);

    % Provincias "frontera": alta variabilidad de estado (0.10 < pi_h < 0.90)
    idx_frontera = find(pi_h_emp > 0.10 & pi_h_emp < 0.90);

    resultados_terr(e).nombre              = escenarios{e};
    resultados_terr(e).p_real_global       = p_real_terr;
    resultados_terr(e).pAceptacion         = mean(aceptaciones);
    resultados_terr(e).pK5_empirica        = mean(K_vals >= 5);
    resultados_terr(e).pK5_teorica         = pK5_teorica;
    resultados_terr(e).K_medio             = mean(K_vals);
    resultados_terr(e).pi_h                = pi_h_emp;
    resultados_terr(e).provincias_frontera = provincias(idx_frontera);

    % --- Salida por escenario ---
    fprintf('\n----------------------------------------------------------------\n');
    fprintf(' ESCENARIO: %s\n', escenarios{e});
    fprintf('----------------------------------------------------------------\n');
    fprintf(' p global real (ponderado)          : %.5f\n', p_real_terr);
    fprintf(' P(aceptacion ILP) empirica          : %.3f\n', resultados_terr(e).pAceptacion);
    fprintf(' P(K>=5) empirica                    : %.3f\n', resultados_terr(e).pK5_empirica);
    fprintf(' P(K>=5) teorica (Poisson-Binomial)   : %.3f\n', pK5_teorica);
    fprintf(' K medio (provincias que superan)    : %.2f / 9\n', resultados_terr(e).K_medio);
    fprintf(' Frecuencia de superacion por provincia (pi_h):\n');
    for h = 1:9
        marca = '';
        if ismember(h, idx_frontera)
            marca = '  <-- FRONTERA';
        end
        fprintf('   %-12s: %.3f%s\n', provincias{h}, pi_h_emp(h), marca);
    end
end

%Tabla resumen comparativa
fprintf('\n\n');
fprintf('================================================================\n');
fprintf('   TABLA RESUMEN - SIMULACION TERRITORIAL CASTILLA Y LEON\n');
fprintf('================================================================\n');
fprintf('%-15s %10s %12s %12s %12s %8s\n', ...
    'Escenario','p_global','P(Acepta)','P(K>=5)emp','P(K>=5)teo','K_medio');
for e = 1:length(escenarios)
    r = resultados_terr(e);
    fprintf('%-15s %10.4f %12.3f %12.3f %12.3f %8.2f\n', ...
        r.nombre, r.p_real_global, r.pAceptacion, r.pK5_empirica, r.pK5_teorica, r.K_medio);
end

% Grafico: distribucion de K (provincias que superan el umbral) por escenario
figure('Name','Simulacion Territorial CyL','NumberTitle','off');
colores_esc = {[0.2 0.5 0.9], [0.85 0.2 0.2], [0.9 0.6 0]};
hold on;
ancho = 0.27;
for e = 1:length(escenarios)
    conteoK = histcounts(K_vals_all(:,e), -0.5:9.5) / Rb_terr;
    offset  = (e - 2) * ancho;
    bar((0:9) + offset, conteoK, ancho, 'FaceColor', colores_esc{e}, ...
        'DisplayName', escenarios{e});
end
xline(4.5, 'k--', 'LineWidth', 1.5, 'Label', 'Umbral mayoria (K>=5)', ...
      'LabelVerticalAlignment','top');
xlabel('K (numero de provincias que superan su umbral)');
ylabel('Probabilidad empirica');
title('Distribucion de K por escenario - Simulacion territorial CyL');
legend('Location','northwest');
grid on;
hold off;