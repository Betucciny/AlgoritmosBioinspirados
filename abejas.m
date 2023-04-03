clear all
format longG
% Numero de fuentes
Nf = 50;
% Numero de variables
Nvar = 3;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = ones(1, Nvar)* -10;
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = ones(1, Nvar)* 10; 
%Numero de iteraciones
Niter = 10000;
%Limite de explotacion de fuente
limite = round(Niter / (2*Nf));

rng('shuffle');
fuentes = crearPob(Li, Ls, Nf, Nvar);
colFO = Nvar + 1;

FO = zeros(Niter, 1);
L = zeros(Niter, 1);
for i = 1:Nf
    FO(i) = funcionObjetivo(fuentes(i,:));
end
best = fuentes(1,:);
bestFO = FO(1);

for p=1:Niter
    asignacion = randperm(Nf);
    for i=1:Nf
        posible = zeros(1,Nvar);
        k = asignacion(i);
        for j=1:Nvar
            phi = -1 + 2 * rand();
            posible(j) = ajustar(fuentes(i,j) + phi*(fuentes(i,j)-fuentes(k,j)), Li(j), Ls(j)); 
        end
        FOposible = funcionObjetivo(posible);
        if FOposible < FO(i)
            fuentes(i,:) = posible;
            FO(i) = FOposible;
            L(i) = 0;
        else
            L(i) = L(i) + 1;
        end
    end
    
    %Asignacion de abejas al azar (todas las abejas con una fuente)
    asignacion = randperm(Nf);
    for i=1:Nf
        posible = zeros(1,Nvar);
        k = asignacion(i);
        for j=1:Nvar
            phi = -1 + 2 * rand();
            posible(j) = ajustar(fuentes(i,j) + phi*(fuentes(i,j)-fuentes(k,j)), Li(j), Ls(j)); 
        end
        FOposible = funcionObjetivo(posible);
        if FOposible < FO(i)
            fuentes(i,:) = posible;
            FO(i) = FOposible;
            L(i) = 0;
        else
            L(i) = L(i) + 1;
        end
    end
    
    for i=1:Nf
        if FO(i) < bestFO
            best = fuentes(i,:);
            bestFO = FO(i);
        end
    end

    for i=1:Nf
        if L(i) > limite
            L(i) = 0;
            fuentes(i,:) = crearPob(Li,Ls, 1, Nvar);
            FO(i) = funcionObjetivo(fuentes(i,:));
        end
    end
    
    p
    best
    bestFO
end

function FO = funcionObjetivo(reales)
    FO = sum(reales.^2);
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