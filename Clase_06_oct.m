%% Ejemplo 1
x1 = [0 0 1 1 1 1]; % tren de impulsos
n1 = [0 1 2 3 4 5];

% ¿Cómo lo representamos? Pista: Sumatoria

%% Ejemplo 2
n2 = -10 : 10;
x21 = n2 == 0; % Muestra unitaria
x22 = n2 >= 3; % Escalón unitario
x23 = n2.*(n2 >= 0); % Rampa unitaria
x24 = 0.7.^(n2); % Exponencial con a = 0.7
x25 = (0.88*exp(1i*pi/8)).^n2; % Exponencial complejo
x26 = 0.99*sin(2*pi*n2/7);

subplot(321), stem(n2,x21); ylabel('Muestra');
subplot(322), stem(n2,x22); ylabel('Escalón');
subplot(323), stem(n2,x23); ylabel('Rampa'); xlabel('Muestras')

subplot(324), stem(n2,x24); ylabel('Exponencial');
subplot(325), stem(n2,abs(x25)); ylabel('Exponencial compleja');
subplot(326), stem(n2,x26); ylabel('Sinusoidal');

%% Ejemplo 3
n3 = -5 : 5;
x3 = abs(n3).*(abs(n3) <= 3);

y31 = x3; % Sistema identidad
y32 = [x3(2:end),0]; % Sistema con retardo unitario
y33 = [0,x3(1:end-1)]; % Sistema con adelanto unitario
y34 = (y32 + y31 + y33)/3; % Sistema de media móvil
y35 = (y32 + y31 + y33); % Sistema de mediana móvil
y36 = sum(tril(ones(numel(x3),1)*x3),2)'; % Acumulador

% diferenciador, modulador, multiplicador, folder
%%
x1 = [0 1 1 0];
x2 = [0 1 1 0];

y = conv(x1,x2);
t = 0:length(y)-1;
stem(t,y);