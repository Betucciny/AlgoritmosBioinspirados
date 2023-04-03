clear all
% Numero de individuos
Nind = 40;
% Numero de variables
Nvar = 4;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0 0 -0.55 -0.55];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [1200 1200 0.55 0.55];
%Numero de generaciones del genetico
Ngen = 20000;
%Factor de cruza
Fc = 0.5;
% Fc = 0.6 + (0.3)*rand();
%Factor de mutacion
Fm = 0.7;


rng('shuffle');
poblacion = crearPob(Li, Ls, Nind, Nvar);
colFO = Nvar + 1;

FO = zeros(Nind, 1);
S = zeros(Nind, 1);
for i = 1:Nind
    FO(i) = funcionObjetivo(poblacion(i,:));
    g = restdes(poblacion(i,:));
    h = restigu(poblacion(i,:));
    S(i) = SVR(g, h);
end

% archivoFO = fopen("bin/binDEBfo(30).txt", "w");
% archivoSR = fopen("bin/binDEBsr(30).txt", "w");
% final = fopen("bin/binDEB(30).txt", "w");

for p = 1:Ngen
    p;
    u = zeros(size(poblacion));
    

    
    % Impresion de resultados

    impr = [poblacion, FO, S];
    sortrows(impr, colFO);
    FOm = impr(1,colFO);
    Sm = impr(1,colFO+1);
    impr(1,1:colFO-1)
    FOm
    Sm
%     fprintf(archivoFO, '%10f ', FOm);
%     fprintf(archivoSR, '%10f ', Sm);

    Fm = 0.5 + (0.3)*rand();
    %Generacion de ruido
    for i = 1:Nind
        indices = randperm(Nind, 3);
        while any(indices==i)
            indices = randperm(Nind, 3);
        end
        u(i, :) = poblacion(indices(1),:) + Fm * (poblacion(indices(3),:) ...
            - poblacion(indices(2),:));
        for elem = 1:Nvar
            u(i,elem) = ajustar(u(i,elem), Li(elem), Ls(elem));
        end
    end

    newPob = zeros(size(poblacion));

    %Creacion de trial y eleccion
    for i = 1:Nind
        jrand = randi(Nvar);
        trial = zeros(1, Nvar);
        for j = 1:Nvar
            r = rand();
            if r < Fc || j==jrand
                trial(1, j) = u(i,j);
            else
                trial(1, j) = poblacion(i,j);
            end
        end

        FOtarget = FO(i);
        Starget = S(i);

        FOtrial = funcionObjetivo(trial);
        Strial = SVR(restdes(trial), restigu(trial));

        %         Reglas de DEB
        if DEB(FOtrial, Strial, FOtarget, Starget)
            newPob(i, :) = trial;
            FO(i) = FOtrial;
            S(i) = Strial;
        else
            newPob(i, :) = poblacion(i, :);
        end
        
    end 
    poblacion = newPob;    
end


% fprintf(final, '%10f %10f %10f %10f %10f %10f', impr(1,:));
% fclose(final);
% fclose(archivoFO);
% fclose(archivoSR);


function FO = funcionObjetivo(p)
    FO = 3*p(1) + 0.000001*p(1).^3 + 2*p(2) + (0.000002/3)*p(2).^3;
end

function s = SVR(g, h)
    s = 0;
    for i = 1:size(g,2)
        s = s + max([0 g(i)]);
    end
    for i = 1:size(h,2)
        s = s + max([0 abs(h(i))-0.0001]);
    end
end

function g = restdes(p)
    g = zeros(1,2);
    g(1) = -p(4) + p(3) - 0.55;
    g(2) = -p(3) + p(4) - 0.55;
end

function h = restigu(p)
    h = zeros(1,3);
    h(1) = 1000*sin(-p(3)-0.25) + 1000*sin(-p(4) - 0.25) + 894.8 - p(1);
    h(2) = 1000*sin(p(3)-0.25) + 1000*sin(p(3) - p(4) - 0.25) + 894.8 - p(2);
    h(3) = 1000*sin(p(4)-0.25) + 1000*sin(p(4) - p(3) - 0.25) + 1294.8;
end

function mejor1= DEB(FO1, SVR1, FO2, SVR2)
FOtrial = FO1;
FOtarget = FO2;
Strial = SVR1;
Starget = SVR2;
if Starget == 0 && Strial == 0
    if FOtrial < FOtarget
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
elseif Starget ~= 0 && Strial ~= 0
    if Starget > Strial
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
else
    if Strial == 0
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
end
end

function pob = crearPob(li, ls, Nind, Nvar)
    pob = zeros(Nind, Nvar);
    for i=1:Nvar
        pob(:,i) = li(i) + (ls(i)- li(i))*rand(Nind, 1);
    end
end

function ajustado = ajustar(valor, li, ls)
    while true
        if valor < li
            valor = 2*li - valor;
        elseif valor > ls
            valor = 2*ls - valor;
        else
            break
        end
    end
    ajustado = valor;
end