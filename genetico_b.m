
clear all
% Numero de individuos
Nind = 20;
% Numero de variables
Nvar = 10;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = ones(1, Nvar)* -10;
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = ones(1, Nvar)* 10;
% Arreglo de tamaño Nvar con las precisiones correspondientes
Press = ones(1, Nvar)* 2;
%Numero de generaciones del genetico
Ngen = 5000;
%Factor de cruza
Fc = 0.85;
%Factor de mutacion
Fm = 0.7;



rng('shuffle');

[poblacion, nbits, NbitsT] = crearPob(Li, Ls, Press, Nind, Nvar);
FO = funcionObjetivo(pob2real(nbits, Li, Ls, poblacion));
poblacion = [poblacion , FO];
nbits;


archivo = fopen("bin extra.txt", "w");



for i= 1:Ngen
    i
    colFO = NbitsT + 1;
    sumaFO = sum(poblacion(:, colFO));
    freqA = zeros(Nind,1);
    Acum = 0;

    %Calculo de frecuencia relativa acumulado
    for j = 1:Nind
        Acum = poblacion(j,colFO)/sumaFO + Acum;
        freqA(j) = Acum; 
    end
    clear Acum
    
    %Seleccion de n padres
    select = zeros(Nind, NbitsT);
    for j = 1:Nind
        r = rand();
        indice = 0;
        for h = 1:Nind
            if freqA(h, 1) > r
                indice = h;
                break
            end
        end
        select(j,:) = poblacion(indice, 1:NbitsT);
    end
    clear freqA indice
    
    %Generacion de hijos
    newPob = zeros(Nind, NbitsT);
    puntoC = round(NbitsT * Fc);
    for j = 1:2:Nind
        newPob(j, 1:puntoC) = select(j, 1:puntoC);
        newPob(j+1, 1:puntoC) = select(j+1, 1:puntoC);
        newPob(j, puntoC+1:NbitsT) = select(j+1, puntoC+1:NbitsT);
        newPob(j+1, puntoC+1:NbitsT) = select(j, puntoC+1:NbitsT);
    end
    clear puntoC
    
    %Mutacion de hijos
    for j = 1:Nind
        for k = 1:NbitsT
            if rand() > Fm
                if newPob(j,k) == 0
                   newPob(j,k) = 1;
                else
                    newPob(j,k) = 0;
                end
            end
        end
    end

    %Seleccion de la siguiente generacion
    FO = funcionObjetivo(pob2real(nbits, Li, Ls, newPob));
    tempnewPob = [newPob , FO];
    temppoblacion = [poblacion ; tempnewPob];
    temppoblacion = sortrows(temppoblacion, colFO, "descend");
    poblacion = temppoblacion(1:Nind, :);

    FO = temppoblacion(1,colFO);
    fprintf(archivo, '%10f ', 1./FO);

end

reales = pob2real(nbits, Li, Ls, poblacion(:,1:NbitsT))
FO = 1./funcionObjetivo(reales)


fprintf(archivo, '\n%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n', reales(1,:));
fprintf(archivo, '%10f ', FO(1));
fclose(archivo);


% Evalucacion de la funcion Objetivo a partir de la matriz de los reales
function FO = funcionObjetivo(reales)
    FO = 1./sum(reales.^2, 2);
end

function real = bits2range(bits, li, ls, nbits)
    real = li+(bits*(ls-li))/(2^nbits - 1);
end


function pobr = pob2real(nbits, li, ls, poblacion)
    nin = size(poblacion, 1);
    nvar = size(nbits, 1);
    pobr = zeros(nin, nvar);
    actual = 1;
    for i=1:nvar
        temp1 = poblacion(:,actual:actual+nbits(i)-1);
        for j = 1:nin
            row = temp1(j,:);
            row_str = bin2dec(num2str(row));
            pobr(j,i) = bits2range(row_str, li(i), ls(i), nbits(i));
        end
        actual = actual+nbits(i);
    end
end

% La poblacion va x1, x2, ...
function [pob, nbits, nbitsT] = crearPob(li, ls, pres, Nind, Nvar)
    nbits = zeros(Nvar,1);
    nbitsT = 0;
    for i=1:Nvar
        nbits(i) = real2nbits(li(i),ls(i),pres(i));
        nbitsT = nbitsT + nbits(i);
    end
    pob = randi([0 1], Nind, nbitsT);
end

function nbits = real2nbits(li,ls,p)
nbits = round(log((ls-li)*10^p)+0.9);
end

