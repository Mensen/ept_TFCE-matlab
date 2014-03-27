
function [] = ERPErrorPlot(x,y,e)

nL = size(y,1); % Specify the number of lines in the plot

if isempty(x)
    x = 1:1:size(y,2);
    x = x(ones(nL,1),:);
end

uL = y+e; %Upper Limit
lL = y-e; %Lower Limit

yP = [lL, fliplr(uL)];
xP = [x,fliplr(x)];

%% Plotting

H.Figure = figure;
set(H.Figure,...
    'Name',             'ERP Comparison'    ,...
    'NumberTitle',      'off'               ,...
    'Units',            'pixels'            ,...
    'Position',         [100 100 900 600]   ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            );

H.Axes = axes;
set(H.Axes,...
    'ColorOrder',       flipud(hsv(nL)),     ...     % set the new ColorOrder to the first three scaled colours of HSV colormap
    'NextPlot',         'replacechildren'   );

H.Line = plot(H.Axes, x',y'); % Create the first line and assign a handle
set(H.Line,...
    'LineWidth',        3                   );

hold on

plot(x',uL',':');
plot(x',lL',':');

H.Patch = patch(xP',yP',1);
set(H.Patch,...
    'EdgeColor',        'none'              ,...
    'FaceColor',        'flat'              ,...
    'FaceVertexCData',  flipud(hsv(nL))     ,...
    'CDataMapping',     'direct'            ,...
    'FaceAlpha',        0.1                 );

end
