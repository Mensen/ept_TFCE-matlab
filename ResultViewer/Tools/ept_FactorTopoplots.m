function ept_FactorTopoplots(Data, e_loc, Results, Factor, Sample)
%% Levels and Statistics Topoplot Tight

% This file is part of the program ept_ResultViewer.
% ept_ResultViewer is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% ept_ResultViewer is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with ept_ResultViewer.  If not, see <http://www.gnu.org/licenses/>.

% Input
% If channel*frequency*time data then the Sample variable should be
% [frequencyBin+size(Data{1},3)*(timeSample-1);

%% Process the Data

if isstruct(Results.TFCE_Obs)
    if Factor == 1
        a       = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
        b       = cellfun(@mean, a, 'UniformOutput', false); 
        Stats   = Results.TFCE_Obs.A(:,Sample)';
    elseif Factor == 2
        Data    = Data';    
        a       = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
        b       = cellfun(@mean, a, 'UniformOutput', false); 
        Stats   = Results.TFCE_Obs.B(:,Sample)';
    else
        fprintf('Only main effects currently available: Factor = 1 or 2.');
        return
    end
else
    a       = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
    b       = cellfun(@mean, a, 'UniformOutput', false);
    b       = b(:);
    Stats =  Results.TFCE_Obs(:,Sample)';
end
%% Calculate the average ERP over multiple conditions
nData = zeros(size(b,1), numel(b{1}));
for i = 1:size(b,2);
    
    nData = nData + cell2mat(b(:,i));
    
end
nData = nData/size(b,2);

%% Calculate plot specific variables

numPlots    = size(nData,1)+1;                     % How many plots in the figure? Each level plus the final statistics plot
axesPos     = (0:1:numPlots)/numPlots;             % Determines the starting x-axis point for each plot (last element redundant)
axesWidth   = 1/numPlots;
mapLimits   = [min(nData(:)), max(nData(:))];      % Calculates the total ERP amplitude range over all groups to adjust the colormaps later

%% Prepare the figure
H.Figure = figure;
set(H.Figure,...
    'Name',             ['S:' num2str(Sample) ' Factor Topoplots']  ,...
    'Units',            'normalized'        ,...    % Normalises the figure size to fit the screen
    'Position',         [0 0 1 0.75]        ,...    % Creates a figure half the screen height
    'NumberTitle',      'off'               ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            ,...
    'DefaultAxesLineStyleOrder',    {'-'}   );

H.MainAxes = axes('Position',[0 0 1 1]);
set(H.MainAxes,...
    'Visible',          'off'               );

% Prepare the individual topoplot axes
for i = 1:numPlots
    
    H.Axes(i) = axes('Position',[axesPos(i) 0.25 axesWidth 0.5]);
    H.Titles(i) = title(H.Axes(i), ['Level ', int2str(i)]);

end

% Adjusts the maps color scheme so topoplot colours are equal among level topoplots
set(H.Axes(1:end-1),...
    'CLim',             mapLimits           );


H.Titles(end) = title('TFCE Stats');
set(H.Titles,...
    'FontSize'   , 16               , ...
    'FontName'   , 'AvantGarde'     ,...
    'FontWeight' , 'bold'           );

%% Plot the actual amplitude topoplots for each group

for i = 1:numPlots-1
   
    set(H.Figure,'CurrentAxes',H.Axes(i))
    ept_Topoplot(nData(i,:), e_loc);
    
end

%% Statistics Topoplot

set(H.Figure,'CurrentAxes',H.Axes(end))
ept_Topoplot(Stats, e_loc);
