clear all
% Numero de individuos
Nind = 30;
% Numero de variables
Nvar = 3;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = ones(1, Nvar)* -5;
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = ones(1, Nvar)* 5;
% Numero de generaciones del genetico
Ngen = 500;
% Inercia
W = 0.1;
% Cognitivo
c1 = 2;
% Social
c2 = 1.2;

rng('shuffle');
colFo = Nvar + 1;

pob = crearPob(Li, Ls, Nind, Nvar)
vel = zeros(size(pob));
FO = zeros([Nind 1]);
for j = 1:Nvar
   vel(:, j) = ones(Nind, 1) * (Li(j) + Ls(j)) / 2;
end
for i = 1:Nind
    FO(i) = funcionObjetivo(pob(1,:));
end
pbest = pob;
pbestFO = FO;

temp = [pob, FO];
sortrows(temp, colFo);
gbest = temp(1,1:Nvar);
gbestFO = temp(1,colFo);

archivo = fopen("pso(29).txt", "w");


for p = 1:Ngen
    newvel = zeros(size(vel));
    newpob = zeros(size(pob));
    newFO = zeros(size(FO));
    

    for i= 1:Nind
        for var = 1:Nvar
            newvel(i, var) = W*vel(i, var) + c1*rand()*(pbest(i,var) - ...
                pob(i,var)) + c2*rand()*(gbest(1,var) - pob(i,var));
            newpob(i, var) = ajustar(pob(i,var) + newvel(i, var), ...
                Li(var), Ls(var));
        end
        newFO(i) = funcionObjetivo(newpob(i,:));

        if newFO(i) < pbestFO(i)
            pbest(i,:) = newpob(i,:);
            pbestFO(i) = newFO(i);
        end
    end

    for i =1:Nind
        if pbestFO(i) < gbestFO
            gbest(1,:) = pbest(i,:);
            gbestFO = pbestFO(i);
        end
    end

    FO = newFO;
    vel = newvel;
    pob = newpob;
    fprintf(archivo, '%10f ', gbestFO);
    gbest
    gbestFO
end

fprintf(archivo, '\n%10f %10f %10f\n', gbest);
fprintf(archivo, '%10f ', gbestFO);
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


