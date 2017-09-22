function my_PSO(fun,fronteras)

%%%%%%%%%%%%%%%%%Preámbulo%%%%%%%%%%%%%%
D = size(fronteras,1);%Dimensiones
M = 100;%Iteraciones
N = 30;%Partículas

C1 = 2.5;%coef. Conf. Enjambre
C2 = 2.5;%coef. Conf. Propia
W  = 0.7;%Coef. Inercia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V  = zeros(N,D);

%%%%%%%%%%%%%Inicialización%%%%%%%%%%%%%%%%%%%%
Aleatorios = rand(N,D);

P  = nan(size(Aleatorios));

for n = 1:N
    P(n,:) = fronteras(:,1)' + (fronteras(:,2)-fronteras(:,1))' .* Aleatorios(n,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(P(:,1),P(:,2),'*')%Prueba

%%%%%%%%%%%%%%%Evaluación de f%%%%%%%%%%%%%%%%%%%
F = nan(N,1); % Función evaluada

for n = 1:N
    F(n) = fun(P(n,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[X,Y] = meshgrid(linspace(fronteras(1,1),fronteras(1,2),50),linspace(fronteras(2,1),fronteras(2,2),50));%Prueba
%surf(X,Y,reshape(fun([X(:),Y(:)]),50,50));shading interp; hold on; colormap winter;
%plot3(P(:,1),P(:,2),F,'ro','MarkerFaceColor','y');

%%%%%%%%%%%%%%%%%%Identificar%%%%%%%%%%%%5
Pgi = P;%Son iguales
[Fbest,g] = min(F);
Pg = P(g,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot3(Pg(1),Pg(2),Fbest,'Or','MarkerFaceColor','r','MarkerSize',10);

%%%%%%%%%%%%%%%%%%%%%%%calculo del movimiento%%%%%%%%%%%%%%%%5
for m = 1:M
    for n = 1:N
        V(n,:) = W*(V(n,:) + C1*rand(1,D) .* (Pg-P(n,:))+...
            C2*rand(1,D) .* (Pgi(n,:)-P(n,:)));
        P(n,:) = P(n,:) + V(n,:);
    end
   
    F(n) = fun(P(n,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Identificar Pgi%%%%%%%%%%%%%%
    [Fbest,g] = min(F);
    Pg = P(g,:);
    
    for n = 1:N
        if fun(Pgi(n,:) > F(n))
            Pgi(n,:) = P(n,:);
        end
    end
    [X,Y] = meshgrid(linspace(fronteras(1,1),fronteras(1,2),50),linspace(fronteras(2,1),fronteras(2,2),50));%Prueba
    surf(X,Y,reshape(fun([X(:),Y(:)]),50,50));shading interp; hold on; colormap winter;
    plot3(P(:,1),P(:,2),F,'ro','MarkerFaceColor','y');
    plot3(Pg(1),Pg(2),Fbest,'Or','MarkerFaceColor','r','MarkerSize',10);
    axis([fronteras(1,:),fronteras(2,:)]);title(sprintf('ite=%d, Fbest = %.4g',m,Fbest));
    hold off;view(2);getframe(gcf);pause(0.05);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end