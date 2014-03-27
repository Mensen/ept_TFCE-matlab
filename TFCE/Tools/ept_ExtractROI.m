function [Output] = ept_ExtractROI(Channels,S1,S2)
%% ept Function to average a signal over listed channels and samples for potential further analysis (e.g. ANCOVAs)

% Load the file
[File, FilePath] = uigetfile('', 'Please Select the TFCE Results File', 'MultiSelect', 'off');
Loaded   = load ([FilePath, '\' File], 'Data', 'Info');
Data     = Loaded.Data;
e_loc    = Loaded.Info.Electrodes.e_loc;

% Specify the channels and sample range
if nargin == 0;
    Input     = inputdlg({'Channel Names', 'Range Start', 'Range End'}, 'Specify Data Range', 1, {'E1,E2,E3', '1', '10'});
    Channels  = regexp(Input{1},',','split');
    S1        = str2num(Input{2});
    S2        = str2num(Input{3});
end

% Find the channels
ChId = zeros(1,numel(e_loc));
for i = 1:numel(Channels)
    ChId   = ChId + strcmp(Channels{i},{e_loc(:).labels});
end
ChId = logical(ChId);

% Get the data...
nData = cellfun(@(x) x(:,ChId,S1:S2), Data, 'UniformOutput', false);
% Average the 2nd and 3rd dimension
mData = cellfun(@(x) mean(x,3), nData, 'UniformOutput', false);
Output = cellfun(@(x) mean(x,2), mData, 'UniformOutput', false);

end



