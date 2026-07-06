function [decision, n_parada, trayectoria] = prps(X_pob, A, B, lr1, lr0)
%   [decision, n_parada, trayectoria] = PRPS(X_pob, A, B, lr1, lr0)
%   Entradas:
%     X_pob       : vector poblacion de 0s y 1s
%     A           : umbral inferior en log-escala (aceptar H0)
%     B           : umbral superior en log-escala (rechazar H0)
%     lr1         : log(p2/p1)       incremento si X_i = 1
%     lr0         : log((1-p2)/(1-p1)) incremento si X_i = 0
%   Salidas:
%     decision    :  1 = Aceptar H0
%                   -1 = Rechazar H0
%                    0 = Indeciso (agotada la poblacion)
%     n_parada    : firmas inspeccionadas hasta la decision
%     trayectoria : vector de log(Lambda_n), longitud n_parada+1

    N_pob       = length(X_pob);
    idx         = randperm(N_pob);        % Orden aleatorio de inspeccion
    log_lam     = 0;
    trayectoria = zeros(N_pob + 1, 1);
    trayectoria(1) = 0;

    for i = 1:N_pob
        xi      = X_pob(idx(i));
        log_lam = log_lam + xi*lr1 + (1-xi)*lr0;
        trayectoria(i+1) = log_lam;
        if log_lam <= A                   % Suficiente evidencia a favor de H0
            decision    =  1;
            n_parada    = i;
            trayectoria = trayectoria(1:i+1);
            return;
        end
        if log_lam >= B                   % Suficiente evidencia contra H0
            decision    = -1;
            n_parada    = i;
            trayectoria = trayectoria(1:i+1);
            return;
        end
    end
    % Si se agota la poblacion sin decision
    decision    =  0;
    n_parada    = N_pob;
    trayectoria = trayectoria(1:N_pob+1);
end
