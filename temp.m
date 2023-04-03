clear clc
clear
% Parámetros del algoritmo
N = 50; % Número de fuentes de alimento
D = 3; % Variables
max_generaciones = 10000; % Número máximo de generaciones
lim_inf = [-5 -5 -5]; % Límite inferior
lim_sup = [5 5 5]; % Límite superior
num_abandonos = round(max_generaciones / (2*N));

% Inicialización de las fuentes de alimento
posiciones = darpob(lim_sup,lim_inf, N, D);

% Evaluación de las fuentes de alimento
valores = zeros(N,1);

for i = 1:N
    valores(i) = funcion_objetivo(posiciones(i,:));
end

% Encontrar la mejor solución global (mejor_posicion)
[mejor_valor, idx] = min(valores);
mejor_posicion = posiciones(idx, :);

% Inicialización del vector de intentos sin mejora
intentos_sin_mejora = zeros(N,1);

% Ciclo de búsqueda
for generacion = 1:max_generaciones

    %Pirmer ciclo
    for i = 1:N
        % Selección aleatoria de otra fuente de alimento
        k = randi([1,N]);

        while k == i
            k = randi([1,N]);
        end

        % Generación de una nueva solución
        phi = -1 + 2*rand(1,D);
        nueva_posicion = posiciones(i,:) + phi.*(posiciones(i,:) - posiciones(k,:));
        nueva_posicion = max(nueva_posicion, lim_inf);
        nueva_posicion = min(nueva_posicion, lim_sup);
        nuevo_valor = funcion_objetivo(nueva_posicion);

        % Actualización de la fuente de alimento si la nueva solución es mejor
        if nuevo_valor < valores(i)
            posiciones(i,:) = nueva_posicion;
            valores(i) = nuevo_valor;
                if nuevo_valor < mejor_valor
                    mejor_posicion = nueva_posicion;
                    mejor_valor = nuevo_valor;
                end
            intentos_sin_mejora(i) = 0; % Reinicia el contador de intentos sin mejora para la fuente i
        else
            intentos_sin_mejora(i) = intentos_sin_mejora(i) + 1; % Incrementa el contador de intentos sin mejora para la fuente i
        end   
   end

    %Segundo ciclo
    for i = 1:N

        % Asigna una abeja a cada fuente de alimento
        k = randi([1,N]);
        while k == i
            k = randi([1,N]);
        end

        % Generación de una nueva solución
        phi = -1 + 2*rand(1,D);
        nueva_posicion = posiciones(i,:) + phi.*(posiciones(i,:) - posiciones(k,:));
        nueva_posicion = max(nueva_posicion, lim_inf);
        nueva_posicion = min(nueva_posicion, lim_sup);
        nuevo_valor = funcion_objetivo(nueva_posicion);

        % Actualización de la fuente de alimento si la nueva solución es mejor
        if nuevo_valor < valores(i)
            posiciones(i,:) = nueva_posicion;
            valores(i) = nuevo_valor;
            if nuevo_valor < mejor_valor
                mejor_posicion = nueva_posicion;
                mejor_valor = nuevo_valor;
            end
            intentos_sin_mejora(i) = 0; % Reinicia el contador de intentos sin mejora para la fuente i
        else
            intentos_sin_mejora(i) = intentos_sin_mejora(i) + 1; % Incrementa el contador de intentos sin mejora para la fuente i
        end
    end
    
    for i = 1:N
        % Intento de abandono de la fuente de alimento
        if intentos_sin_mejora(i) >= num_abandonos
            % Genera una nueva solución aleatoria dentro de los límites
            posiciones(i,:) = lim_inf + (lim_sup-lim_inf).*rand(1,D);
            valores(i) = funcion_objetivo(posiciones(i,:));
            % Actualiza la mejor solución si es necesario
            if valores(i) < mejor_valor
                mejor_posicion = posiciones(i,:);
                mejor_valor = valores(i);
            end
            intentos_sin_mejora(i) = 0; % Reinicia el contador de intentos sin mejora para la fuente i
        end
    end

    disp(['Generación ', num2str(generacion), ': mejor valor = ', num2str(mejor_valor)]);
end

% Función objetivo (sphere)
function valor = funcion_objetivo(x)
    valor = sum(x.^2);
end


function [pob] = darpob(Nls,Nli, pobl, Nvariables)

    pob=zeros(pobl,Nvariables);
    for i = 1:pobl
       for j = 1:Nvariables
           pob(i,j) = (Nls(j)-Nli(j))*rand() + Nli(j);
      end
    end
end