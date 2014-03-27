function ept_WaveletChannelOverview(Data, Results, Ch, Thresh)
%% Plots channel Wavelets for each level with the TFCE overview

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

%% Process the Input

if isstruct(Results.TFCE_Obs)
%     if Factor == 1
%         a       = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
%         b       = cellfun(@mean, a, 'UniformOutput', false); 
%         Stats   = Results.TFCE_Obs.A(:,Sample)';
%     elseif Factor == 2
%         Data    = Data';    
%         a       = cellfun(@(x) x(:,:,Sample), Data, 'UniformOutput', false);
%         b       = cellfun(@mean, a, 'UniformOutput', false); 
%         Stats   = Results.TFCE_Obs.B(:,Sample)';
%     else
        error('Currently only available for direct T-Test comparisons...')
%     end
else
    a       = cellfun(@(x) x(:,Ch,:,:), Data, 'UniformOutput', false);
    b       = cellfun(@mean, a, 'UniformOutput', false);
    nData   = cellfun(@squeeze, b, 'UniformOutput', false);
    nData   = nData(:);
    Stats   = squeeze(Results.TFCE_Obs(Ch,:,:));
    P_Values= squeeze(Results.P_Values(Ch,:,:));
    nData   = cellfun(@(x) mag2db(x), nData, 'UniformOutput', false); %Transform the data to dB (just remove if already like that).
end

%% Calculate plot specific variables
numPlots    = size(nData,1);                     % How many plots in the figure? Each level plus the final statistics plot
axesPos     = (0:1:numPlots)/numPlots;             % Determines the starting x-axis point for each plot (last element redundant)
axesWidth   = 1/numPlots;
mapLimits = [min(min([nData{:}])), max(max([nData{:}]))]; % Calculates the total ERP amplitude range over all groups to adjust the colormaps later

%% Prepare the figure
H.Figure = figure;
set(H.Figure,...
    'Name',             ['Wavelet Summary for Channel ' num2str(Ch)]  ,...
    'Units',            'normalized'        ,...    % Normalises the figure size to fit the screen
    'Position',         [0 0 0.5 0.5] ,...    % Creates a figure half the screen height
    'NumberTitle',      'off'               ,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            ,...
    'DefaultAxesLineStyleOrder',    {'-'}   );

% Prepare the level topoplot axes
for i = 1:numPlots
    
    H.Axes(i) = axes('Position'                 ,...
                [axesPos(i)+0.05                ,...
                numPlots/(numPlots+1)+0.02      ,... % initial height proportional to number of plots
                axesWidth-0.08                  ,...
                1/(numPlots+1)-0.08             ...  % total height dependent to number of plots
                ]);
    H.Titles(i) = title(H.Axes(i), ['Level ', int2str(i)]);

end

% Main Stats Axes
H.Axes(end+1) = axes('Position'             ,...
                [0+0.05                      ,...
                0+0.1                      ,...
                1-0.075                      ,...
                numPlots/(numPlots+1)-0.2  ...  % total height dependent to number of plots
                ]);

% Title and labels
H.Titles(end+1) = title('TFCE Stats');
H.XLabel = xlabel('Time (sample)');
H.YLabel = ylabel('Frequency (bin)');
set(H.Titles,...
    'FontSize'   , 16               , ...
    'FontWeight' , 'bold'           );
set([H.Titles,H.XLabel, H.YLabel],...
    'FontName'   , 'AvantGarde'     );
            
% Adjusts the maps color scheme so topoplot colours are equal among level topoplots
set(H.Axes(1:end-1),...
    'CLim',             mapLimits           );

%% Plot the actual amplitude topoplots for each grou

for i = 1:numPlots
   
    set(H.Figure,'CurrentAxes',H.Axes(i))
    ept_WaveletPlot(nData{i});
    
end

%% Statistics Topoplot

set(H.Figure,'CurrentAxes',H.Axes(end))
ept_WaveletPlot(Stats, 'PValues', P_Values, 'LevelList', Thresh);




