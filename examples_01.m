%% Example -1

plot(rand(1000,1),'*k'), ylim([-0.5 1.5]);

randab = @(n,a,b) a + (b - a)*rand(n,1);

hold on,
plot(randab(1000,-3,5),'*r'), ylim([-3.5 5.5]);

%% Example 0

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

% X = S + 2*randn(size(t)); % awgn ?
X = awgn(S,0.01,'measured');

stem(t,X,'r'), hold on,
stem(t,S,'k'),

figure,

Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Example 1 - Calculating and plotting a DTFT

% multiple signals a^n * u[n] a<1

%format compact;
a       = 0.88 * exp(2i*pi/5);
nn      = 0:40;
xn      = a.^nn;

%plot3(nn,real(xn),imag(xn))


% Obtain the DTFT
[X,W]   = DTFT(xn,2^7);

% Get components
W_normalised    = W/2/pi;
X_magnitude     = abs(X);
X_phase         = 180/pi*angle(X);

% Plotting process
hf          = figure('Color','w');
hf.Name     = 'Example 1 - Calculating and plotting a DTFT';
hf.Color    = 'White';
hf.Units    = 'Normalized';
hf.Position = [0.3 0.2 0.43 0.6];

% First plot
subplot(211), h1 = plot(W_normalised, X_magnitude,'k','LineWidth',2);
title('Magnitude response','Interpreter','LaTeX');
ylabel('$$| H (\omega) |$$','Interpreter','LaTeX','FontSize',18);
xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

subplot(212), h2 = plot(W_normalised, X_phase,'k','LineWidth',2);
title('Phase response','Interpreter','LaTeX');
ylabel('$$\phi_{H} (\omega)$$ [deg]','Interpreter','LaTeX','FontSize',18);
xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);
h2.Parent.TickLabelInterpreter = 'LaTeX';
h2.Parent.FontSize  = 18;
h2.Parent.LineWidth = 1.5;

%% Example 2 - Calculating and plotting a DTFT

t       = linspace(0,40,200); % s ?
T       = diff(t(1:2));
L       = numel(t);
nn      = 0 : L - 1;
xn      = 0.89 * sin(2*pi*T*nn/4) + 0.32 * cos(2*pi*T*nn/3);

% Obtain the DTFT
[X,W]   = DTFT(xn, 2^9, true);

%% Example 3 - Calculating and plotting a DTFT (Dirichlet kernel)

nn      = 0:50;
xn      = (nn <= 25);

% Obtain the DTFT
[X,W]   = DTFT(xn, 2^9, true);

%% HW: write a asinc(w,L) function and compare its results
 % modify DTFT to shiftDTFT(x,n0,N) to find the DFT of x[n-n0]
 % what does freqz do?
%% Example 4 - Flipping the frequency axis

a       = 0.88 * exp(2i*pi/5);
nn      = 0:40;
xn      = a.^nn;

% Obtain the DTFT
[X,W]   = DTFT(xn, 2^9);

% check
[Xflipped, Wflipped] = flipDTFT(X, W);

% plot
hf = figure('Color','w');
hf.Units    = 'Normalized';
hf.Position = [0.3 0.2 0.43 0.6];

h1 = plot(W,abs(X),'k','LineWidth',2); hold on,
h2 = plot(Wflipped,abs(Xflipped),'--r','LineWidth',2);

ylabel('$$| H (\omega) |$$','Interpreter','LaTeX','FontSize',18);
xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);

h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;
legend([h1,h2],{'DTFT','flipped DTFT'},'Interpreter','LaTeX','FontSize',18)

%% Example 5 - Flipping the frequency axis

nn      = -20 : 20; % s ?
L       = 10;
xn      = (L - nn).*(nn >= 0 & nn <= L) + (L + nn).*(nn < 0 & nn > -L);

% Obtain the DTFT
[X,W]   = DTFT(xn, 2^9);

% check
[Xflipped, Wflipped] = flipDTFT(X, W);

% plot
hf = figure('Color','w');
hf.Units    = 'Normalized';
hf.Position = [0.3 0.2 0.43 0.6];

h1 = plot(W,abs(X),'k','LineWidth',2); hold on,
h2 = plot(Wflipped,abs(Xflipped),'--r','LineWidth',2);

ylabel('$$| H (\omega) |$$','Interpreter','LaTeX','FontSize',18);
xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);

h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;
legend([h1,h2],{'DTFT','flipped DTFT'},'Interpreter','LaTeX','FontSize',18)

%% Example 5 - Convergence
clear all;

% Infinitely long exponential
N   = 500;
nn  = 0 : N - 1;
a   = 0.977;

% Create the signal functional
xnt = a.^(nn);                  % true
xnw = @(L) a.^(nn).*(nn < L);   % windowed

% True DTFT
[Xt,Wt]   = DTFT(xnt,2^9);

% Plot with different values of L
hf = figure('Color','w');
hf.Name     = 'Signals';
hf.Units    = 'Normalized';
hf.Position = [0.3 0.1 0.43 0.8];

L = 32; xn(1,:)  = xnw(L); 
subplot(221), h1 = plot(nn,xn(1,:),'k','LineWidth',2); 
hold on, plot(nn,xnt,'r--','LineWidth',2); 
line([L L],[0 1],'Color','r','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(L)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

L = 64; xn(2,:)  = xnw(L);
subplot(222), h1 = plot(nn,xn(2,:),'k','LineWidth',2);
hold on, plot(nn,xnt,'r--','LineWidth',2); 
line([L L],[0 1],'Color','r','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(L)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

L = 128; xn(3,:)  = xnw(L);
subplot(223), h1 = plot(nn,xn(3,:),'k','LineWidth',2);
hold on, plot(nn,xnt,'r--','LineWidth',2); 
line([L L],[0 1],'Color','r','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(L)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

L = 256; xn(4,:)  = xnw(L);
subplot(224), h1 = plot(nn,xn(4,:),'k','LineWidth',2);
hold on, plot(nn,xnt,'r--','LineWidth',2); 
line([L L],[0 1],'Color','r','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(L)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

% DTFTs
hf = figure('Color','w');
hf.Name     = 'Signals';
hf.Units    = 'Normalized';
hf.Position = [0.3 0.1 0.43 0.8];

[X,W]   = DTFT(xn(1,:),2^9);
subplot(221), h1 = semilogy(W/2/pi,abs(X),'k','LineWidth',2);
hold on, h2 = semilogy(Wt/2/pi,abs(Xt),'r--','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(32)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

[X,W]   = DTFT(xn(2,:),2^9);
subplot(222), h1 = semilogy(W/2/pi,abs(X),'k','LineWidth',2);
hold on, h2 = semilogy(Wt/2/pi,abs(Xt),'r--','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(64)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

[X,W]   = DTFT(xn(3,:),2^9);
subplot(223), h1 = semilogy(W/2/pi,abs(X),'k','LineWidth',2);
hold on, h2 = semilogy(Wt/2/pi,abs(Xt),'r--','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(128)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;

[X,W]   = DTFT(xn(4,:),2^9);
subplot(224), h1 = semilogy(W/2/pi,abs(X),'k','LineWidth',2);
hold on, h2 = semilogy(Wt/2/pi,abs(Xt),'r--','LineWidth',2);
ylabel('$$x[n]$$','Interpreter','LaTeX','FontSize',18);
xlabel('$$ n $$','Interpreter','LaTeX','FontSize',18);
title(['$$ L $$=',num2str(256)],'Interpreter','LaTeX','FontSize',18);
h1.Parent.TickLabelInterpreter = 'LaTeX';
h1.Parent.FontSize  = 18;
h1.Parent.LineWidth = 1.5;
