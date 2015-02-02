function [output] = ept_ExtractROI(channels, sample_start, sample_end)
%% ept Function to average a signal over listed channels and samples for potential further analysis (e.g. ANCOVAs)

% Load the file
[file_name, file_path] = uigetfile('', 'Please Select the TFCE Results File', 'MultiSelect', 'off');
Loaded   = load (fullfile(file_path, file_name), 'Data', 'Info');
Data     = Loaded.Data;
e_loc    = Loaded.Info.Electrodes.e_loc;

% Specify the channels and sample range
if nargin == 0;
    Input     = inputdlg({'Channel Names', 'Range Start', 'Range End'}, 'Specify Data Range', 1, {'E1,E2,E3', '1', '10'});
    channels  = regexp(Input{1},',','split');
    sample_start        = str2num(Input{2});
    sample_end        = str2num(Input{3});
end

% Find the channels
if isa('channels', 'cell')
    ChId = zeros(1,numel(e_loc));
    for i = 1:numel(channels)
        ChId   = ChId + strcmp(channels{i},{e_loc(:).labels});
    end
    ChId = logical(ChId);
else
    ChId = channels;
end

% Get the data...
nData = cellfun(@(x) x(:, ChId, sample_start:sample_end), Data,...
    'UniformOutput', false);

% Average the 2nd and 3rd dimension
mData = cellfun(@(x) mean(x,3), nData, 'UniformOutput', false);
output = cellfun(@(x) mean(x,2), mData, 'UniformOutput', false);

end



