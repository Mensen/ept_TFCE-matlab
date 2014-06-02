function ept_ERPWaves(Info, Data, Results, Channels)

factor = 1;
bL = 200;

% Channels = {'E78','E79', 'E85', 'E86', 'E87', 'E92', 'E93'}; % Right Parietal
% Channels = {'E24', 'E20', 'E27', 'E28','E29', 'E34', 'E35'}; % Left Central
% Channels = {'E14','E8', 'E9', 'E1', 'E2'}; % Right Frontal
% Channels = {'E74', 'E75', 'E76', 'E81', 'E82', 'E83', 'E88', 'E89'}; % Right Occipital

if nargin == 3;
    Input     = inputdlg({'Channel Names', 'Range Start', 'Range End'}, 'Specify Data Range', 1, {'E1,E2,E3', '1', '10'});
    Channels  = regexp(Input{1},',','split');
%     S1        = str2num(Input{2});
%     S2        = str2num(Input{3});
end

e_loc    = Info.Electrodes.e_loc;

% Find the channels
if isa(Channels, 'double')
    ChId = Channels;
else
    ChId = zeros(1,numel(e_loc));
    for i = 1:numel(Channels)
        ChId   = ChId + strcmp(Channels{i},{e_loc(:).labels});
    end
    ChId = logical(ChId);
end

if ~isfield(Info.Parameters, 'GroupSizes')
    Info.Parameters.GroupSizes = cell2mat(cellfun(@(x) size(x, 1), Data, 'UniformOutput', false));
end

if factor == 1;
    for i = 1:size(Data,1)
        
        m  = squeeze(mean(cell2mat(cellfun(@(x) x(:,ChId,:), Data(i,:), 'UniformOutput', false)')));
        y(i,:) = mean(m);
        s  = squeeze(std(cell2mat(cellfun(@(x) x(:,ChId,:), Data(i,:), 'UniformOutput', false)')))/sqrt(Info.Parameters.GroupSizes(1));
        e(i,:) = mean(s);
        
    end
elseif factor == 2;
    for i = 1:size(Data,2)
        % Calculate SE for within group using Cousineau method
        ParticipantScore = cell2mat(cellfun(@(x) x(:,ChId,:), Data(:,i), 'UniformOutput', false));
        ParticipantMean  = mean(cell2mat(cellfun(@(x) x(:,ChId,:), Data(:,i), 'UniformOutput', false)),2);
        ParticipantMeanRep = repmat(ParticipantMean, [1,size(ParticipantScore,2),1]);
        ConditionDifferences = ParticipantScore-ParticipantMeanRep;
        s = squeeze(mean(ConditionDifferences));
        e(i,:) = mean(s);
        m   = squeeze(mean(cell2mat(cellfun(@(x) x(:,ChId,:), Data(:,i), 'UniformOutput', false))));
        y(i,:) = mean(m);        
    end
end

x = (1:1:(Info.Parameters.nSamples))*(1000/Info.Parameters.rSample)-bL;

nL = size(y,1); % Specify the number of lines in the plot
x = x(ones(nL,1),:);

uL = y+e; %Upper Limit
lL = y-e; %Lower Limit

yP = [lL, fliplr(uL)];
xP = [x,fliplr(x)];

%% Actual ERP Plotting

H.Figure = figure;
set(H.Figure,...
    'Name',             'ERP Comparison'    ,...
    'NumberTitle',      'off'               ,...
    'Units',            'pixels'            ,...
    'Position',         [100 100 900 500]   ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            ,...
    'DefaultAxesColorOrder',        flipud(hsv(nL))       ,...
    'DefaultAxesLineStyleOrder',    {'-',':', '--', '-.'} );
H.Axes = axes;
set(H.Axes,...
    'FontName',         'Helvetica'         ,...
    'XLim',             [x(1,1),x(1,end)]   ,...
    'NextPlot',         'replacechildren'   );
    
hold on;

H.Line = plot(H.Axes, x',y'); % Create the first line and assign a handle
set(H.Line,...
    'LineWidth',        3                   );

Y_Lim = get(H.Axes, 'YLim');
line([0,0],Y_Lim, 'Color', [0 0 0], 'LineStyle', '--')

plot(x',uL',':');
plot(x',lL',':');

H.Patch = patch(xP',yP',1);
set(H.Patch,...
    'EdgeColor',        'none'              ,...
    'FaceColor',        'flat'              ,...
    'FaceVertexCData',  flipud(hsv(nL))     ,...
    'CDataMapping',     'direct'            ,...
    'FaceAlpha',        0.15                );

c = arrayfun(@(col) strcat(Channels{1,col}, ', '), 1:size(Channels,2),'Unif', false);

hTitle  = title (['ERP Comparison for ', c{:}]);
hXLabel = xlabel('Time (ms)'                     );
hYLabel = ylabel('ERP Strength (MicroVolts)'     );

for i = 1:nL
    Legend{i} = ['Level ' int2str(i)];
end
H.Legend = legend(H.Line, Legend);
        
