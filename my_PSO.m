function my_PSO(fun,fronteras)

%Preambulo
D=size(fronteras,1);
M=100;
N=30;
V=zeros(N,D);
C1=2.5;
C2=2.5;
W=0.7;

%Inicializacion
Aleatorios=rand(N,D);
P=nan(size(Aleatorios));
for n=1:N
    P(n,:)=fronteras(:,1)' + (fronteras(:,2)-fronteras(:,1))'.*Aleatorios(n,:);
end

%plot(P(:,1),P(:,2),'o')
%>>my_PSO(o,[-10,10;-5,5])

%Evaluacion de funcion
F=nan(N,1);
for n=1:N
    F(n)=fun(P(n,:));
end

%Prueba 2D
[X,Y]= meshgrid(linspace(fronteras(1,1),fronteras(1,2),50),linspace(fronteras(2,1),fronteras(2,2),50));
surf(X,Y,reshape(fun([X(:),Y(:)]),50,50));
shading interp; hold on; colormap winter;
plot3(P(:,1),P(:,2),F,'ro','MarkerFaceColor','y');

%>>my_PSO(@(X) X(:,1).^2 + X(:,2).^2, [-5 5, -10 10]);

Pgi=P; %Son iguales comenzando 
[Fbest,g]=min(F);
Pg=P(g,:);
plot3(Pg(1),Pg(2),Fbest,'or','MarkerFaceColor','r','Markersize',10);
hold on; getframe(gcf);

%Pgi = El historial de la mejor posicion de la particula

%Calculo de P^(t+1)
for m=1:M
    for n=1:N
        V(n,:)=W*(V(n,:)+C1*rand(1,D).*(Pg-P(n,:))+...
            C2*rand(1,D).*(Pgi(n,:)-P(n,:)));
        P(n,:)=P(n,:)+V(n,:);
    end

    F=nan(N,1);
for n=1:N
    F(n)=fun(P(n,:));
end

%Identificar Pg
[Fbest,g]=min(F);
Pg=P(g,:);

%Identificar Pgi
for n=1:N
    if fun(Pgi(n,:)) > F(n)
        Pgi(n,:)=P(n,:);
    end
end

[X,Y]= meshgrid(linspace(fronteras(1,1),fronteras(1,2),50),linspace(fronteras(2,1),fronteras(2,2),50));
surf(X,Y,reshape(fun([X(:),Y(:)]),50,50));
shading interp; hold on; colormap winter;
plot3(P(:,1),P(:,2),F,'ro','MarkerFaceColor','y');
plot3(Pg(1),Pg(2),Fbest,'or','MarkerFaceColor','r','Markersize',10);
hold off; 
axis([fronteras(1,:),fronteras(2,:)]);
getframe(gcf);
    
end;
