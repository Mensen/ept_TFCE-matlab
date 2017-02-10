function H = ept_IndividualTopoplots(Data, e_loc, Results, Factor, Sample)
% Ept function to plot each individuals topography for a given sample...

% Input
% If channel*frequency*time data then the Sample variable should be
% [frequencyBin+size(Data{1},3)*(timeSample-1)];

if isstruct(Results.TFCE_Obs)
    if Factor == 1
        iData   = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
        b       = cellfun(@mean, iData, 'UniformOutput', false); 
        Stats   = Results.TFCE_Obs.A(:,Sample)';
    elseif Factor == 2
        Data    = Data';    
        iData   = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
        b       = cellfun(@mean, iData, 'UniformOutput', false); 
        Stats   = Results.TFCE_Obs.B(:,Sample)';
    else
        error('Only main effects currently available: Factor = A or B.')
    end
else
    iData  = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
    b       = cellfun(@mean, iData, 'UniformOutput', false);
    Stats =  Results.TFCE_Obs(:,Sample)';
end

% Calculate the average ERP over multiple conditions (just cell2mat for single condition TTest)
gData = zeros(size(b,1), numel(b{1}));
% will have to make multiple condition means for iData too...
for i = 1:size(b,2);
    
    gData = gData + cell2mat(b(:,i));
    
end
gData = gData/size(b,2);


% Calculate plot specific variables
% Find max number of participants in each group and add 3 (TFCE Stats will be twice as large...

group_sizes = cellfun(@(x) size(x,1), Data);
numPlots    = max(group_sizes(:))+3;

axesPos     = (0:1:numPlots)/numPlots;             % Determines the starting x-axis point for each plot (last element redundant)
axesWidth   = 1/numPlots;
mapLimits   = [min(gData(:)), max(gData(:))];      % Calculates the total ERP amplitude range over all groups to adjust the colormaps later

% Prepare the figure
H.Figure = figure;
set(H.Figure,...
    'Name',             ['S:' num2str(Sample) ' Factor Topoplots']  ,...
    'Units',            'normalized'        ,...    % Normalises the figure size to fit the screen
    'Position',         [0 0 1 0.6]         ,...    % Creates a figure half the screen height
    'NumberTitle',      'off'               ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            ,...
    'DefaultAxesLineStyleOrder',    {'-'}   );

% Prepare the individual topoplot axes
for j = 1:size(iData,1) % loop for rows (levels)
    for i = 1:numPlots-2 % loop for plots (individuals + group)

        H.Axes(i,j) = axes(...
            'parent', H.Figure,...
            'Position', [axesPos(i) (j-1)/size(iData,1) axesWidth 1/size(iData,1)] ,...
            'xtick', [] ,...
            'ytick', [] );
%         H.Titles(i,j) = title(H.Axes(i), ['Level ', int2str(i)]);

    end
end

% Draw the larger TFCE Stats axes
H.TFCEAxes = axes('Position', [axesPos(end-2), 0, axesWidth*2, 1]);

% Adjusts the maps color scheme so topoplot colours are equal among level topoplots
set(H.Axes,'CLim', mapLimits);

% Plot the actual topoplots of each individual
for j = 1:size(iData,1) % loop for rows (levels)
    for i = 1:size(iData{j},1) % loop for plots (individuals + group)
        % draw each topoplot
        ept_Topoplot(iData{j}(i,:), e_loc,...
            'Axes', H.Axes(i,j), ...
            'PlotChannels', 0, ...
            'LineStyle', 'none');
    end
end 

% Plot group topoplots
for i = 1:size(gData,1)
    ept_Topoplot(gData(i,:), e_loc,...
        'Axes', H.Axes(end, i), ...
        'PlotChannels', 0, ...
        'LineStyle', 'none');
end

% Plot the larger TFCE topoplot
ept_Topoplot(Stats, e_loc,...
    'Axes', H.TFCEAxes, ...
    'LineStyle', 'none');
