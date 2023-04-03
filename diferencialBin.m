clear all
% Numero de individuos
Nind = 40;
% Numero de variables
Nvar = 10;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = ones(1, Nvar)* -10;
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = ones(1, Nvar)* 10;
%Numero de generaciones del genetico
Ngen = 200;
%Factor de cruza
Fc = 0.5;
%Factor de mutacion
Fm = 0.7;

rng('shuffle');
poblacion = crearPob(Li, Ls, Nind, Nvar)

archivo = fopen("dif(1).txt", "w");

for p = 1:Ngen
    p;
    colFO = Nvar + 1;
    u = zeros(size(poblacion));
    FO = funcionObjetivo(poblacion);

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
        FOtrial = funcionObjetivo(trial);
        if FOtrial < FO(i)
            newPob(i, :) = trial;
        else
            newPob(i, :) = poblacion(i, :);
        end
    end
    poblacion = newPob;

    impr = [poblacion, funcionObjetivo(poblacion)];
    sortrows(impr, colFO);
    FOm = impr(1,colFO);
    impr(1,1:colFO-1)
    FOm
    fprintf(archivo, '%10f ', FOm);

end

fprintf(archivo, '\n%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n', impr(1,1:colFO-1));
fprintf(archivo, '%10f ',impr(1,colFO));
fclose(archivo);

function FO = funcionObjetivo(reales)
    FO = sum(reales.^2, 2);
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