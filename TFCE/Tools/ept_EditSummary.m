function ept_EditSummary(Channel)
%% Review ERP summary files for participant outliers

[File, DataPath] = uigetfile('', 'Please Select the EEG Summary File', 'MultiSelect', 'off');
Data     = load ([DataPath, '\' File], 'Summary', 'DataFile');
Summary  = Data.Summary;
Labels   = Data.DataFile;

if nargin == 0;
    Channel  = 121;
end

H.Figure = figure;
set(H.Figure,...
    'Name',             ['Explore Individual ERPs']  ,...
    'Units',            'normalized'        ,...    % Normalises the figure size to fit the screen
    'Position',         [0 0 1 1]           ,...    % Creates a figure the full screen height
    'NumberTitle',      'off'               ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            );

H.CurrentAxes = axes;
set(H.CurrentAxes,...
    'Position',         [0.1 0.1 1 0.9]   );

H.ERP_Plot = plot(1:size(Summary,3), squeeze(Summary(:,Channel,:)));

H.Legend    = legend(H.ERP_Plot, Labels);
set(H.Legend,...
    'Visible',          'off'             );
% set(H.Legend,...
%     'Location',         'EastOutside'     );

plottools(H.Figure, 'on', 'plotbrowser')

end