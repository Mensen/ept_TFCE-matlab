%% Find the difference between two summary files and save the outpit
% Primarily useful for then calculating a three-way interaction despite
% using only two factors (as 1 factor is reduced to a difference)

function ept_AverageSummaries

clear
clc

[DataFiles, DataPath] = uigetfile('*.mat', ['Select the summary files to average'], 'MultiSelect', 'on');

for i = 1:length(DataFiles)
    
    LoadData = load(fullfile(DataPath, DataFiles{i}), 'Summary');
    Data{i}         = LoadData.Summary;
end

% Check for equal size data files


% Calculate the mean of the datasets...
dim = ndims(Data{1});          %# Get the number of dimensions for your arrays
M = cat(dim+1,Data{:});        %# Convert to a (dim+1)-dimensional matrix
Summary = mean(M,dim+1);       %# Get the mean across arrays

saveName = [DataFiles{1,1}(1:(end-4)), 'Avg'];

prompt = {'Desired Savename?'};
def = {saveName};
noGr = inputdlg(prompt, 'User Input', 1, def);

saveName=noGr{1};

save ([saveName, '.mat'], 'Summary', 'DataFiles', '-mat')

end