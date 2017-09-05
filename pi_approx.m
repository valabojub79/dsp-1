% Numero de puntos
N   = 10000;

% Generar puntos
dardos = rand(N,2);

% Approximacion
PI = 4*sum(sqrt(sum(dardos.^2,2)) <= 1)/N;