function [fig] = adefig (opt, fig, axi)
%ADEFIG     Adequate Figure.
% ADEFIG is a simple code intended to format a common ugly figure.
% 
% [fig] = ADEFIG (opt, fig, axi)
%
% Inputs:   opt - option control with two possible options (currently).
%           opt : 'begin'   - pre-format a non-existing figure yet.
%           opt : 'end'     - format fig existing figure with (at least)
%                             its axes (axi).
%           fig - figure object. If fig does not exist as an object, 
%                 ADEFIG will take the current figure.
%           axi - axes object. If axi does not exist as an object, 
%                 ADEFIG will take the current axes.
%
% Example:
%
% x = linspace(-pi/2,pi/2,500);
% hf = ADEFIG('begin');             % Pre-format
% ha = plot(x,sin(10*x),'k','LineWidth',1.5); ylim([-1.5 1.5]);
% ADEFIG('end', hf, ha.Parent);     % Format
%
% See also [put here other codes from Jorge's toolbox].
%
% Copyright 2017 - Jorge Mario Cruz-Duarte (mrcrois@gmail.com).

% If opt is not given, show help
if nargin < 1, 
    help adefig;
    return;
end

% Properties
Interpreter             = 'LaTeX';
FontSize                = 18;
LineWidth               = 1.5;


switch opt,
    case 'begin',
        % Create the figure
        fig = figure;
        
        % Set figure properties
        fig.Name        = 'This figure has formatted using adefig.m';
        fig.Color       = 'White';
        fig.Units       = 'Normalized';
        fig.Position    = [0.3 0.2 0.4 0.6];
        fig.PaperUnits  = 'centimeters';
        fig.PaperSize   = [60 60];

    case 'end',
        if nargin < 2, fig = get(gcf); end
        if nargin < 3, axi = get(gca); end
        
        axi.TickLabelInterpreter    = Interpreter;
        axi.FontSize                = FontSize;
        axi.LineWidth               = 1;%LineWidth;
        axi.Box                     = 'on';
        
        if numel(axi.XLabel.String) == 0, axi.XLabel.String = 'x'; end
        axi.XLabel.Interpreter      = Interpreter;
        axi.XLabel.FontSize         = FontSize;
        axi.XMinorTick              = 'on';
        
        if numel(axi.YLabel.String) == 0, axi.YLabel.String = 'y'; end
        axi.YMinorTick              = 'on';
        axi.YLabel.Interpreter      = Interpreter;
        axi.YLabel.FontSize         = FontSize;
        
        if axi.View(1) ~= 0,
            if numel(axi.ZLabel.String) == 0, axi.ZLabel.String = 'z'; end
            axi.ZMinorTick          = 'on';
            axi.ZLabel.Interpreter      = Interpreter;
            axi.ZLabel.FontSize         = FontSize;
        end
            
    otherwise,
        error('I don''t know what I should do when you say: %s! :''(',opt);
end