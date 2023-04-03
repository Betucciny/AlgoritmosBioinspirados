clear all
% Numero de individuos
Nind = 20;
% Numero de variables
Nvar = 10;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = ones(1, Nvar)* -10;
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = ones(1, Nvar)* 10;
%Numero de generaciones del genetico
Ngen = 5000;
%Factor de cruza
Fc = 0.85;
%Pocision de cruza
Pc = 0.85;
%Factor de mutacion
Fm = 0.7;

rng('shuffle');
poblacion = crearPob(Li, Ls, Nind, Nvar)
FO = funcionObjetivo(poblacion);
poblacion = [poblacion, FO];


archivo = fopen("reales (30).txt", "w");

for i = 1:Ngen
    i
    colFO = Nvar + 1;
    sumaFO = sum(poblacion(:, colFO));
    freqA = zeros(Nind,1);
    Acum = 0;

    %Calculo de frecuencia relativa acumulado
    for j = 1:Nind
        Acum = poblacion(j,colFO)/sumaFO + Acum;
        freqA(j,1) = Acum;
    end
    clear Acum

     %Seleccion de n padres
    select = zeros(Nind, Nvar);
    for j = 1:Nind
        r = rand();
        for h = 1:Nind
            if freqA(h) > r
                indice = h;
                break
            end
        end
        select(j,:) = poblacion(indice, 1:Nvar);
    end
    clear freqA indice


    %Generacion de hijos
    newPob = ones(Nind, Nvar);
    puntoC = round(Nvar * Pc);
    for j = 1:2:Nind
        newPob(j, 1:puntoC) = select(j, 1:puntoC);
        newPob(j+1, 1:puntoC) = select(j+1, 1:puntoC);
        newPob(j, puntoC+1:Nvar) = Fc .*select(j+1, puntoC+1:Nvar) + (1-Fc).* select(j, puntoC+1:Nvar);
        newPob(j+1, puntoC+1:Nvar) =  Fc .*select(j, puntoC+1:Nvar) + (1-Fc).* select(j+1, puntoC+1:Nvar);
    end
    clear puntoC
    

    %Mutacion de hijos
    for j = 1:Nind
        for k = 1:Nvar
            if rand() > Fm
                li = Li(k);
                ls = Ls(k);
                newPob(j,k) = ajustar(li + (ls- li)*rand(), li, ls);
            end
        end
    end


    %Seleccion de la siguiente generacion

    FO = funcionObjetivo(newPob);
    tempnewPob = [newPob , FO];
    temppoblacion = [poblacion ; tempnewPob];
    temppoblacion = sortrows(temppoblacion, colFO, "descend");
    poblacion = temppoblacion(1:Nind, :);
    poblacion;

    fprintf(archivo, '%10f ', 1./temppoblacion(1, colFO));

end

poblacion(1,1:colFO-1)
1./poblacion(1,colFO)

fprintf(archivo, '\n%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n', poblacion(1,1:colFO-1));
fprintf(archivo, '%10f ', 1./poblacion(1,colFO));
fclose(archivo);


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


function FO = funcionObjetivo(reales)
    FO = 1./sum(reales.^2, 2);
end

function pob = crearPob(li, ls, Nind, Nvar)
    pob = zeros(Nind, Nvar);
    for i=1:Nvar
        pob(:,i) = li(i) + (ls(i)- li(i))*rand(Nind, 1);
    end
end