clear all
% Ireaciones
g_max = 100000;
% Harmony numbers
NH = 20;
% Variable numbers
Nvar = 3;
% Inferior limits
Li = ones(1, Nvar)* -10;
% Upper limits
Ls = ones(1, Nvar)* 10;
% Acceptance rate
rac = 0.9;
% Pitch adjust rate
rpa = 0.3;
% Intelligent acceptance
ria = 0.5;
% Factor para cambiar bw
a = 1;

harmonies = crearPob(Li, Ls, NH, Nvar);
FO = zeros([1 NH]);


% Calculo de las primeras funciones objetivo
for i = 1:NH
    FO(i) = funcionObjetivo(harmonies(i,:));
end


% Encontremos el peor de nuestras harmonias
windex = indicePeor(FO, NH);
bindex = indiceMejor(FO, NH);

for g = 1:g_max
    bw = (Ls - Li)/g.^a;
    newH = zeros([1 Nvar]);
    for v = 1:Nvar
        if rand() < rac
            index = randi(NH);
            if rand() < rpa
                newH(v) = harmonies(index, v) + bw(v) * (-1 + 2*rand());
            else
                if rand() < ria
                    newH(v) = harmonies(bindex, v);
                else
                    newH(v) = harmonies(index, v);
                end
            end
        else
            newH(v) = crearPob(Li(v), Ls(v), 1, 1);
        end
    end
    newFO = funcionObjetivo(newH);

    if newFO < FO(windex)
        harmonies(windex, :) = newH;
        FO(windex) = newFO;
        windex = indicePeor(FO, NH);
        bindex = indiceMejor(FO, NH);
    end
    harmonies(bindex, :);
    FO(bindex);
end

harmonies(bindex, :)
FO(bindex)


function windex = indicePeor(FO, NH)
windex = 1;
for i = 1:NH
    if FO(i) > FO(windex)
        windex = i;
    end
end
end

function bindex = indiceMejor(FO, NH)
bindex = 1;
for i = 1:NH
    if FO(i) < FO(bindex)
        bindex = i;
    end
end
end


function FO = funcionObjetivo(reales)
FO = sum(reales.^2, "all");
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