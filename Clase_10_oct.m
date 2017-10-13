%% Example 1
ns = [0 1 2 3 4 5 6 7 8 9];
xs = [1 1 1 0 0 0 0 0 1 1];

% Sample of a signal
subplot(311), stem(ns,xs);

% Signal with five periods
xr = repmat(xs,1,5);
subplot(312), stem(xr);

% Fourier series coefficients
Xf = fft(xs);
subplot(313), stem(real(Xf));

%% Example 2
ns = [0 1 2 3 4 5 6 7 8 9];
xs = [1 1 1 0 0 0 0 0 1 1];

% Lowpass filter
H = [1 1 1 0 0 0 0 0 1 1];

% Get the FFT of X
Xf = fft(xs);

% Filter x using H
Yf = Xf.*H;

% Obtain x from Xf.*H
y = ifft(Yf);

% Plot solution
stem(ns,xs,'b'); hold on; stem(ns,y,'r');

%% Example 3
dt = 0.001;
t = 0:dt:2*2*pi;
x = cos(t);
subplot(211), plot(t,x);

% Obtain FFT
f = fftshift(fft(x));
N = length(f);
n = -(N-1)/2:(N-1)/2;
w = 2*pi*n/N/dt;
f = f ./ N;
subplot(212), stem(w,real(f));
axis([-10 10 -1 1]);

%% Example 4
dt = 0.001;
t = 0:dt:2.3*2*pi;
x = cos(t);
figure, subplot(211), plot(t,x);

% Obtain FFT
f = fftshift(fft(x));
N = length(f);
n = -(N-1)/2:(N-1)/2;
w = 2*pi*n/N/dt;
f = f ./ N;
subplot(212), stem(w,real(f));
axis([-10 10 -1 1]);

%% Example 5
dt = 0.001;
t = 0:dt:10;
Nx = length(t);
x = [ones(1,floor(N/2)) zeros(1,ceil(N/2))];

% Perform the convolution x*x
y = conv(x,x);
Ny = length(y);
ty = 0:dt:dt*(Ny-1);
subplot(211), plot(ty,y);

% Perform the convolution in frequency
f = fft(x,Ny);
g = f.*f;
z = ifft(g,Ny);
tz = ty;
subplot(212), plot(tz,z);

%% Example 6
dt = 0.001;
tx = 0:dt:2*2*pi;
Nx = length(tx);
x1 = [ones(1,floor(Nx/2)) zeros(1,ceil(Nx/2))];
x2 = cos(tx);

% Time convolution
y = conv(x1,x2);
Ny = length(y);
ty = 0:dt:dt*(Ny-1);
subplot(211), plot(ty,y);

% Freq multiplication
f1 = fft(x1,Ny);
f2 = fft(x2,Ny);
g = f1.*f2;
z = ifft(g);
tz = ty;
subplot(212), plot(tz,z);

%% Example 7
len = 100;
n = 0:1:len-1;
x = zeros(1,len);
x([0:1:4]+1) = 1;

% Plot it
subplot(311), stem(n,x);

% Plot its fft
G = fft(x);
w = linspace(0,2*pi,len);
subplot(312), stem(w,G);

% Changing the number of samples
len = 20;
n = 0:1:len-1;
x = zeros(1,len);
x([0:1:4]+1) = 1;
G = fft(x);
w = linspace(0,2*pi,len);
subplot(313), stem(w,G);

%% Example 8, time shifting
len = 100;
n = 0:1:len-1;
w = linspace(0,2*pi,len);
x = (n <= 4+1);
G = fft(x);

% Time shifting
y = [zeros(1,10) x(1,1:90)]; % n0 = 10
H = fft(y);
F = exp(-1i*w*10).*G;

% Plot them
subplot(211), stem(w,H);
subplot(212), stem(w,F);

%% Example 9
len = 100;
n = 0:1:len-1;
w = linspace(0,2*pi,len);
x = (n <= 4+1);
G = fft(x);

% Frequency shiting
y = exp(-1i*pi/3*n).*x;
H = fft(y);

% Plot them
subplot(211), stem(w,G);
subplot(212), stem(w,H);

%% Example 10
wc = pi/4;
len = 100;
w = linspace(0,2*pi,len);
n = 0:1:len-1;

G = zeros(1,len);

% find the index of the cutoff frequency
wlow = max(find(w<wc));
whigh = min(find(w>(2*pi-wc)));

% set the freq-response to 1 in the proper regions
G(1:wlow) = 1;
G(whigh:end) = 1;

% Plot it
subplot(211), stem(w,G);

% Calculate the IR of the filter
x = ifft(G);

% Multiply by (-1)^n -> h_hp [n] = (-1)^n h_lp [n]
y = x.*((-1).^n);
F = fft(y);
subplot(212), stem(w,real(F));

%% Example 11
t = 0:0.0001:1;
x = sin(3*pi*t) + 2*cos(50*pi*t) + sin(100*pi*t) + sin(77*pi*t);
c = sqrt(2)*cos(2*pi*1000*t);

% Modulation
y = x.*c;

subplot(311), plot(t,x);
subplot(312), plot(t,y);

% Multiplication
r = y.*c;% + 1*randn(size(y));

% Filtering
[B,A] = butter(3,0.01);
d = filter(B,A,r);

% Plot all
subplot(313), plot(t,x,'b',t,d,'r');

%% Example 12
t = 0:0.0001:1;
x = sin(3*pi*t) + 2*cos(50*pi*t) + sin(100*pi*t) + sin(77*pi*t);
c = sqrt(2)*cos(2*pi*1000*t);

% Modulation
y = x.*c;

subplot(311), plot(t,x);
subplot(312), plot(t,y);

% Multiplication
c2 = sqrt(2)*cos(2*pi*1000*t + pi/6);
r = y.*c2;% + 1*randn(size(y));

% Filtering
[B,A] = butter(3,0.01);
d = filter(B,A,r);

% Plot all
subplot(313), plot(t,x,'b',t,d,'r');

%% Example 13
t = 0:0.0001:1;
x = sin(2*pi*t) + 2*cos(5*pi*t);
c = sqrt(2)*cos(2*pi*1000*t);
A = 5;

% Modulation
y = (A+x).*c;

% Demodulation
r = sqrt(0.5*(y.^2 + imag(hilbert(y)).^2));

% Plot
subplot(211), plot(t,r-A,'r',t,x,'b');

% increasing the mdulation index m > 1
A = 1;
y = (A+x).*c;
r = sqrt(0.5*(y.^2 + imag(hilbert(y)).^2));
subplot(212), plot(t,r-A,'r',t,x,'b');

%% Example 14
t = linspace(0,1,500);
x = sin(2*pi*t) + 2*cos(5*pi*t) + rand(size(t));
w = gausswin(50);
y = conv(x,w);

subplot(211), plot(t,x);
subplot(212), plot(t,y(1:end-numel(w)+1));


%% Example 15
I = (rgb2gray(imread('image.png')));
figure, imshow(I);

Ix = imnoise(I,'gaussian');
figure, imshow(Ix);

Iy = im2double(I);

K = @(sig,x,y) exp(-(x.^2+y.^2)/2/sig^2);
[dx,dy] = meshgrid([-2:2]);
weight = K(0.5,dx,dy)/sum(sum(K(0.5,dx,dy)));

Iy=conv2(Iy,weight,'same');
figure, imshow(mat2gray(Iy));


%% Example 16
t = linspace(0,1,500);
x = (t>=0.5) + 0.2*rand(size(t));

N = 10; sig = sqrt(2);
n = linspace(-4,4,N);
w = exp(-n.^2/2/sig^2)/sqrt(2*pi)/sig;
y = conv(x,w);

dw = -n.*w/sig^2;
dy = conv(x,dw);

subplot(311), plot(t,x);
subplot(312), plot(t,y(1:end-numel(w)+1));
subplot(313), plot(t,dy(1:end-numel(dw)+1));


%% Example 16
I = (rgb2gray(imread('image.png')));
figure, imshow(I);

Ix = imnoise(I,'gaussian');
figure, imshow(Ix);

Iy = im2double(I);

K = @(sig,x,y) exp(-(x.^2+y.^2)/2/sig^2);
[dx,dy] = meshgrid([linspace(-5,5,3)]);

% Primero
sigma = 0.25;
weight = K(sigma,dx,dy)/sum(sum(K(sigma,dx,dy)));

dwx = -dx.*weight/sigma^2;
dwy = -dy.*weight/sigma^2;

Iyx=conv2(Iy,dwx,'same');
Iyy=conv2(Iy,dwy,'same');

dIy = sqrt(Iyx.^2 + Iyy.^2);
dIy_norm = (dIy - min(min(dIy)))/( max(max(dIy)) -  min(min(dIy)));

figure, imshow(mat2gray(dIy_norm));
title(sprintf('\\sigma = %.2f',sigma));

% Segundo
sigma = 1.5;
weight = K(sigma,dx,dy)/sum(sum(K(sigma,dx,dy)));

dwx = -dx.*weight/sigma^2;
dwy = -dy.*weight/sigma^2;

Iyx=conv2(Iy,dwx,'same');
Iyy=conv2(Iy,dwy,'same');

dIy = sqrt(Iyx.^2 + Iyy.^2);
dIy_norm = (dIy - min(min(dIy)))/( max(max(dIy)) -  min(min(dIy)));

figure, imshow(mat2gray(dIy_norm));
title(sprintf('\\sigma = %.2f',sigma));

% Tercero
sigma = 2;
weight = K(sigma,dx,dy)/sum(sum(K(sigma,dx,dy)));

dwx = -dx.*weight/sigma^2;
dwy = -dy.*weight/sigma^2;

Iyx=conv2(Iy,dwx,'same');
Iyy=conv2(Iy,dwy,'same');

dIy = sqrt(Iyx.^2 + Iyy.^2);
dIy_norm = (dIy - min(min(dIy)))/( max(max(dIy)) -  min(min(dIy)));

figure, imshow(mat2gray(dIy_norm));
title(sprintf('\\sigma = %.2f',sigma));


%% addition

dIy_norm = (dIy - min(min(dIy)))/( max(max(dIy)) -  min(min(dIy)));
dIy_umb = dIy_norm > 0.12;

figure, imshow(mat2gray(dIy_umb));
title(sprintf('\\sigma = %.2f',sigma));

%% Example 17   
s = 1;

[x,y] = meshgrid([-s 0 s]);

G = @(x,y,sigma) exp(-(x.^2 + y.^2)/2/sigma^2);

wG = round(G(x,y,s)*3);

dwGx = -x.*wG;
dwGy = -y.*wG;

%% Example 18 : 1 to 

% 1. Leer sonido
[f_,fs] = audioread('sonidito.wav');
Ts = 1/fs;
samples_ = numel(f_);
time_ = 0 : Ts : Ts*(samples_-1);

% 2. Tomar una muestra
cliptime = 10;
time = time_(time_ <= cliptime);
s = f_(time_ <= cliptime);
samples = numel(s);

% 3. Reproduce
pOrig = audioplayer(s,fs);
%pOrig.play;

% 4. Muestra la señal y su espectro
figure; 
subplot(211),
stem(time, s); xlabel('time [s]'), ylabel('amplitude')

ds = fs / samples;
w = (-(samples/2):(samples/2)-1)*ds;
y = fft(s, samples) / samples; 
y2 = fftshift(y);
subplot(212), plot(w,abs(y2));
xlabel('frequency [Hz]'), ylabel('amplitude')

% 5. Agregar ruido
s_n  = awgn(s,10,'measured');
pOrig = audioplayer(s_n,fs);
%pOrig.play;

% 6. Muestra la señal ruidosa y su espectro
figure; 
subplot(211),
stem(time, s_n); xlabel('time [s]'), ylabel('amplitude')

ds = fs / samples;
w = (-(samples/2):(samples/2)-1)*ds;
y = fft(s_n, samples) / samples; 
y2 = fftshift(y);
subplot(212), plot(w,abs(y2));
xlabel('frequency [Hz]'), ylabel('amplitude')

