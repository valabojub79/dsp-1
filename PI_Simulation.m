%Formación del cuarto de circulo
t = 0:1/999:pi/2;
x = cos(t);
y = sin(t);
hold on;
plot(x,y,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de puntos
N   = 10000;

% Generar puntos
dardos = rand(N,2);

% Approximacion
PI = 4*sum(sqrt(sum(dardos.^2,2)) <= 1)/N;

for i=1:1:N 
    %%Obtención de las posociones (x,y) para simular los puntos
    plot(dardos(1*i),dardos(2*i),'*');
    
    %%Muestra de pi en la pantalla
    set(gcf,'Name', [ 'Pi = ',num2str(PI)]);
    
    %$Tiempo de muestreo
    getframe(gcf);
end
