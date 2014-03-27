% Creates a valid summary file for TFCE analysis from individual .dat files

function ept_CreateSummary(File_ID)

if nargin == 0
    File_ID = ''; % Use this to specify a particular part of the file name for easy searching
end

[DataFile, DataPath] = uigetfile(['*' File_ID '*'], 'Please Multi-Select the Files to Merge', 'MultiSelect', 'on');

disp(['You selected ' num2str(length(DataFile)) ' files for this summary... '])

for i = 1:length(DataFile);

    data = importdata ([DataPath '\' DataFile{i}], ' ', 0);
    disp(['Importing file (' num2str(i) '), ' DataFile{i} '... '])
    Summary(i,:,:) = data';
    
end

saveName = DataFile{1,1}(1:(end-4));

prompt = {'Desired Savename?'};
def = {saveName};
noGr = inputdlg(prompt, 'User Input', 1, def);

saveName=noGr{1};

save (['Summary_' saveName '.mat'], 'Summary', 'DataFile')

disp(['Summary file was saved as Summary_' saveName '.mat in the data directory']);
% plot(Summary(:,:,125)')

end 