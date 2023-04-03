clear all
format longG
% Numero de individuos
Nind = 80;
% Numero de variables
Nvar = 4;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0 0 -0.55 -0.55];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [1200 1200 0.55 0.55];
% Numero de generaciones del genetico
Ngen = 100000;
% Inercia
W = 0.7;
% Cognitivo
c1 = 0.8;
% Social
c2 = 1.9;

rng('shuffle');
colFo = Nvar + 1;

pob = crearPob(Li, Ls, Nind, Nvar)
vel = zeros(size(pob));

for j = 1:Nvar
   vel(:, j) = ones(Nind, 1) * (Li(j) + Ls(j)) / 2;
end

FO = zeros(Nind, 1);
S = zeros(Nind, 1);
for i = 1:Nind
    FO(i) = funcionObjetivo(pob(i,:));
    g = restdes(pob(i,:));
    h = restigu(pob(i,:));
    S(i) = SVR(g, h);
end

pbest = pob;
pbestFO = FO;
pbestSVR = S;

gbest = pob(1,:);
gbestFO = FO(1);
gbestSVR = S(1);

for i =1:Nind
    if DEB(pbestFO(i), pbestSVR(i), gbestFO, gbestSVR)
        gbest(1,:) = pbest(i,:);
        gbestFO = pbestFO(i);
        gbestSVR = pbestSVR(i);
    end
end


% archivo = fopen("dif(1).txt", "w");



for p = 1:Ngen
    newvel = zeros(size(vel));
    newpob = zeros(size(pob));
    newFO = zeros(size(FO));
    newSVR = zeros(size(S));

    for i= 1:Nind
        for var = 1:Nvar
            newvel(i, var) = W*vel(i, var) + c1*rand()*(pbest(i,var) - ...
                pob(i,var)) + c2*rand()*(gbest(1,var) - pob(i,var));
            newpob(i, var) = ajustar(pob(i,var) + newvel(i, var), ...
                Li(var), Ls(var));
        end

        newFO(i) = funcionObjetivo(pob(i,:));
        g = restdes(pob(i,:));
        h = restigu(pob(i,:));
        newSVR(i) = SVR(g, h);

        if DEB(newFO(i), newSVR(i), pbestFO(i), pbestSVR(i))
            pbest(i,:) = newpob(i,:);
            pbestFO(i) = newFO(i);
            pbestSVR(i) = newSVR(i);
        end

    end

    for i =1:Nind
        if DEB(pbestFO(i), pbestSVR(i), gbestFO, gbestSVR)
            gbest(1,:) = pbest(i,:);
            gbestFO = pbestFO(i);
            gbestSVR = pbestSVR(i);
        end
    end

    FO = newFO;
    S = newSVR;
    vel = newvel;
    pob = newpob;

%     disp(pob(1:5,:))
%     disp(pbestFO(1:5))
%     disp(pbestSVR(1:5))


    gbest
    gbestFO
    gbestSVR
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

function FO = funcionObjetivo(p)
FO = 3*p(1) + 0.000001*p(1).^3 + 2*p(2) + (0.000002/3)*p(2).^3;
end

function s = SVR(g, h)
s = 0;
for i = 1:size(g,2)
    gtemp = g(i);
    s = s + max([0 gtemp]);
end
for i = 1:size(h,2)
    if abs(h(i)) < 0.0001
        htemp = 0;
    else
        htemp = abs(h(i));
    end
    s = s + max([0 htemp]);
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


