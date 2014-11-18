function H = ept_WaveletPlot(data, varargin)
%% - Set defaults - %%
ContourWidth        = 3;
LevelList           = [0.1,0.05,0.01]; % Plot contours at each of these points
PValues             = [];
DBValues            = 0;

NewFigure           = 0;        % Try to put inside current figure

%% - Process Parameters - %%
if nargin > 1
  if (round(nargin/2) == nargin/2)
    error('Even number of input arguments??')
  end
  for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~ischar(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    
    switch Param
        case 'levellist'
            LevelList           = Value; 
        case 'contourwidth'
            ContourWidth        = Value;
        case 'threshold'
            SigThreshold        = Value;
        case 'newfigure'
            NewFigure           = Value;
        case 'pvalues'
            PValues             = Value;     
        otherwise
            display (['Unknown parameter setting: ' Param])
    end
  end
end

%% - Basic Figure - %%

if isempty(get(0,'children')) || NewFigure == 1
    H.Figure = figure;
    set(H.Figure,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            );

%     set(H.Figure,'Menubar','none');

    H.CurrentAxes = axes('Position',[0.10 0.10 0.8 0.8]);
    
else
    H.CurrentAxes = gca;
end

% Prepare the axes
set(H.CurrentAxes,...
    'XLim',             [1, size(data,2)]         ,...
    'YLim',             [1, size(data,1)]         ,...
    'NextPlot',         'add'               );

%% Plot Surface

% H.Surf = surf(H.CurrentAxes, data);
% set(H.Surf,...
%     'FaceColor',        'interp',           ...
%     'EdgeColor',        'none'              ...
%     );

H.Surf = contourf(H.CurrentAxes, data, ...
    'EdgeColor', 'none');
    

%% Plot the contours of significance

% Adjust the contour lines to account for the minimum and maximum difference in values
if ~isempty(PValues)
    
    [~,H.Contour] = contour(H.CurrentAxes, PValues);
    set(H.Contour,...
        'EdgeColor',        'none',              ...
        'LineWidth',        ContourWidth,        ...
        'Color',            'w',                 ...
        'LevelList',        LevelList,           ...
        'HitTest',          'off'                ...
        );
    
end

% uicontrol('Style','pushbutton',...
%   'String','Zoom In',...
%   'Position',[20 20 60 20],...
%   'Callback','if camva <= 1;return;else;camva(camva-1);end');

