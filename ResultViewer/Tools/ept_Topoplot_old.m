%% New Topoplot Function (ept)

% % Load the results file
% sample  = 92;
% a       = cellfun(@(x) x(:,:,sample), Data, 'UniformOutput', false);
% b       = cellfun(@mean, a, 'UniformOutput', false);
% % New dataset for that particular sample
% nData   = (cell2mat(b(:,1))+ cell2mat(b(:,2)))/2;
% e_loc   = Info.Electrodes.e_loc;
% 
% %% Select a group
% 
% V = nData(1,:);

function H = ept_Topoplot(V, e_loc)

%% Use the e_loc to project points to a 2D surface

Th = pi/180*(cell2mat({e_loc.theta}));      % Calculate theta values from x,y,z e_loc
Rd = cell2mat({e_loc.radius});              % Calculate radian values from x,y,z e_loc

x = Rd.*cos(Th);                            % Calculate 2D projected X
y = Rd.*sin(Th);                            % Calculate 2D projected Y

% Squeeze the coordinates into a -0.5 to 0.5 box
intrad = min(1.0,max(abs(Rd))); intrad = max(intrad,0.5); squeezefac = 0.5/intrad;

x = x*squeezefac; y = y*squeezefac;

%% Create the plotting mesh

gridscale = 100;                            % Creates a 100x100 mesh to map onto 
Xq = linspace(-0.5,0.5,gridscale);
XYq = Xq(ones(gridscale,1),:);

%% Create the interpolation function
x=x(:); y=y(:); V=V(:);                     % Ensure data is in column format
F = TriScatteredInterp(x,y,V, 'natural');
Zi = F(XYq', XYq);

%% Actual Plot

% Prepare the figure

% Check if there is a figure currently opened
if isempty(get(0,'children'))
    H.Figure = figure;
    H.CurrentAxes = axes('Position',[0 0 1 1]);
else
    H.Figure = gcf;
    H.CurrentAxes = gca;
end
% set(H.Figure,...
%     'Color',            'w'                 ,...
%     'Renderer',         'openGL'            );
% 
% % Prepare the axes
% H.MainAxes = axes('Position',[0 0 1 1]);
% set(H.MainAxes,...
%     'XLim',             [-0.5, 0.5]         ,...
%     'YLim',             [-0.5, 0.5]         ,...
%     'NextPlot',         'add'               ,...
%     'Visible',          'off'               );

% Plot the contour map
nContours = 15;                                         % Define number of contour lines to be plotted
levelStep = (max(V)-min(V))/nContours;                  % Translate the number of contours in to steps

[~,H.Contour] = contourf(H.CurrentAxes, XYq,XYq',Zi);
set(H.Contour,...
    'EdgeColor',        'none'              ,...
    'LineWidth',        0.3                 ,...
    'Color',            'k'                 ,...
    'LevelStep',        levelStep           );

% Plot the surface interpolation

% unsh = (gridscale+1)/gridscale; % un-shrink the effects of 'interp' SHADING
% 
% H.Surface = surface(XYq*unsh,XYq'*unsh,zeros(size(Zi)),Zi);
% set(H.Surface,...
%     'EdgeColor',        'none'              ,...
%     'FaceColor',        'interp'            );


% Plot the Head, Ears, and Nose (thanks EEGLAB!)

angle=0:1:360;
datax=(cos(angle*pi/180))/3;
datay=(sin(angle*pi/180))/3; 

H.Head(1) = plot(H.CurrentAxes, datax, datay);

sf = 0.333/0.5; %Scaling factor for the headsize

% Nose...
base  = 0.4954;
basex = 0.0900;                 % nose width
tip   = 0.5750; 
tiphw = 0.02;                   % nose tip half width
tipr  = 0.005;                  % nose tip rounding

H.Head(2) = plot(H.CurrentAxes,...
         [basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf);            % plot nose

% Ears...
q = .04; % ear lengthening
EarX  = [.497-.005  .510        .518        .5299       .5419       .54         .547        .532        .510    .489-.005]; % rmax = 0.5
EarY  = [q+.0555    q+.0775     q+.0783     q+.0746     q+.0555     -.0055      -.0932      -.1313      -.1384  -.1199];


H.Head(3) = plot(H.CurrentAxes,...
                EarX*sf,EarY*sf);% plot left ear
H.Head(4) = plot(H.CurrentAxes,...
                -EarX*sf,EarY*sf);   % plot right ear
  
set(H.Head,...
    'Color',            'k'                 ,...
    'LineWidth',        3                   );

% Adjustments
axis square
axis off
colormap(jet)


