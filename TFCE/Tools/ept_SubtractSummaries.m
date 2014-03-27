%% Find the difference between two summary files and save the outpit
% Primarily useful for then calculating a three-way interaction despite
% using only two factors (as 1 factor is reduced to a difference)

function ept_SubtractSummaries

clear
clc

[DataFiles, DataPath] = uigetfile('*.mat', ['Select the two summary files to subtract'], 'MultiSelect', 'on');

for i = 1:length(DataFiles)
    
    LoadData        = load (strcat(DataPath, DataFiles{i}));
    Data{i}         = LoadData.Summary;
end

Summary = Data{1}-Data{2};

saveName = [DataFiles{1,1}(1:(end-4)), 'D'];

prompt = {'Desired Savename?'};
def = {saveName};
noGr = inputdlg(prompt, 'User Input', 1, def);

saveName=noGr{1};

save ([saveName, '.mat'], 'Summary', 'DataFiles', '-mat')

end