function [H, W] = DTFT (h, N, SHOW)
%DFTF calculate DTFT at N equally spaced frequencies
%   usage:  H = DTFT( h, N)
%       h:  finite-length input vector, whos length is L
%       N:  number of frequencies for evaluation over [-pi,pi]
%               ===> constraint: N >= L
%    SHOW:  if it is TRUE, magnitude and phase are show 
%           	default SHOW = FALSE)
%
%       H: DTFT values (complex)
%       W: (2nd output) vector of freqs where DTFT is computed
%
if nargin < 3, SHOW = false; end
N   = fix(N);
L   = length(h);    h = h(:);   % <-- for vectors ONLY !!!

if (N < L)
    error(['DTFT: the number of data samples cannot exceed ',...
        'the number of frequency samples']);
end

W           = (2*pi/N) * (0:(N - 1))';
mid         = ceil(N/2) + 1;
W(mid:N)    = W(mid:N) - 2*pi;  % <-- move [pi,2*pi) to [-pi,0)

W           = fftshift(W);
H           = fftshift(fft(h,N));   % <-- move negative frequency components

if SHOW == true,
    % Get components
    W_normalised    = W/2/pi;
    H_magnitude     = abs(H);
    H_phase         = 180/pi*angle(H);

    % Plotting process
    hf          = figure('Color','w');
    hf.Name     = 'Example 1 - Calculating and plotting a DTFT';
    hf.Color    = 'White';
    hf.Units    = 'Normalized';
    hf.Position = [0.3 0.2 0.43 0.6];

    % First plot
    subplot(211), h1 = plot(W_normalised, H_magnitude,'k','LineWidth',2);
    title('Magnitude response','Interpreter','LaTeX');
    ylabel('$$| H (\omega) |$$','Interpreter','LaTeX','FontSize',18);
    xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);
    h1.Parent.TickLabelInterpreter = 'LaTeX';
    h1.Parent.FontSize  = 18;
    h1.Parent.LineWidth = 1.5;

    subplot(212), h2 = plot(W_normalised, H_phase,'k','LineWidth',2);
    title('Phase response','Interpreter','LaTeX');
    ylabel('$$\phi_{H} (\omega)$$ [deg]','Interpreter','LaTeX','FontSize',18);
    xlabel('Normalised frequency, $$ \omega/{2\pi} $$','Interpreter','LaTeX','FontSize',18);
    h2.Parent.TickLabelInterpreter = 'LaTeX';
    h2.Parent.FontSize  = 18;   
    h2.Parent.LineWidth = 1.5;
end