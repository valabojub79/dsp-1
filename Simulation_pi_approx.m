clc
clear all
close all

n  = 0;       
N  = 1000;    %Numero de tiros    
DC = 0;       %Numero de dardos dentro del circulo
DF = 0;       %Numero de dardos fuera del circulo

D = rand(N,2);  %Genera la matrix con la posicion de los dardos

%Se dibuja un cuarto de circulo de radio 1
ang = 0:0.0001:pi/2;
x = 1*cos(ang);
y = 1*sin(ang);
plot(x,y,'b','linewidth',1)
hold on, grid on 


while (n < N)
    
d(n+1) = sqrt((D(n+1,1)^2)+(D(n+1,2)^2)); %Calcula la distancia de 
                                          %cada punto respecto al origen
if (d(n+1) <= 1)
    DC = DC+1; colour = 'r';              
else
    DF = DF+1; colour = 'k';
end
n = n+1;

set(gcf,'name',sprintf('pi = %.4f',4*DC/n)); %Muestra en la figura el valor
                                             %aproximado de pi segun el
                                             %numero de iteracion en tiempo
                                             %real
plot(D(n,1),D(n,2),'*','color',colour);
xlabel('Distancia en X')
ylabel('Distancia en Y')
title('Dardos')
axis([0 1 0 1])
getframe(gcf);
end

disp('Numero de dardos dentro del circulo ')
disp(DC)

disp('Numero de dardos fuera del circulo')
disp(DF)
    


