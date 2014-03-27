%% EEG Permutation Test ANOVA (ept) using Threshold-Free Cluster Enhancement 
% Copyright(C) 2012  Armand Mensen (14.12.2010)

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% [Description]
% This tests initially computes a F-Value for each channel and sample
% These values are then enhanced or reduced depending on the size of
% the F-value and the size of the cluster the channel belongs to at 
% various thresholds...
%     
% TFCE is essentially the sum of the clusters to the power of 0.5 multiplied
% by the height of the current threshold squared
% 
% [Input]
% Participant data should be organised into two levels
%     e.g. Dataset 1 is controls while Dataset 2 is patients
% Electrode locations file created using eeglab is required to 
% calculate channel neighbours
% 
% Analysis Types
% i = independent sample T-Test
% d = dependent (paired) sample T-Test
% c = Pearson's correlation (r)
% m = mixed ANOVA
% r = repeated measures ANOVA
%
%
% This Summary file should be a "Participants x Channels x Samples" variable
% - Participants should be in the same order for each factor is repeated
% - Channels must be in the same order as the corresponding electrodes file
% - Samples should be in chronological order
% - Values in the variable correspond to the mV averaged ERP data
% 
% [Output]
% - Info Structure
% - Data Structure
% - Results Structure
%

%% Revision History
%
% 18.02.2012
% Added checks for number of channels in data and locations file
%
% 23.12.2012
% - New permutation strategy for ept_mixedANOVA (2-step model)
% 
% 21.12.2012
% Rewrote the data input windows for more than two conditions
% Have to rewrite the permutation for general amount of levels per factor!
%
% 15.10.2012
% Output is a Info/Results Structure
% Compatible with three groups (still only 2 conditions per factor)
%
% 10.03.2011
% - Now a function that can be called from the command line
%
% 24.05.2011
% - Added random stream dependence on the clock
% 
% 12.07.2011
% - Adapted to loading the two summary files to be compared
% - Uses new modified version of the channel neighbours algorithm

% 10.04.2012
% - Rearranged Data Input and Selection Help

% 09.05.2012
% - Changed default E_H parameters to [0.66, 1] since we use F-values not T


 function [] = ept_TFCE_ANOVA

%% Set Defaults
 
 E_H        = [0.666, 1]; % default parameters of E and H (H is 1 since we are using F-tests; the square of T)
 nPerm      = 2500; % default number of permutations
 rSample    = 250; % default assumed sampling rate (not used in calculation but only for result viewer)
 saveName   = ['ept_Results_', date, '.mat'];
 type       = 'm'; % 'r' = repeated measures; 'm' = mixed ANOVA
 nGroups 	= 2;
 plots      = 0; % change to '1' to show significance plots after calculation
  
% Set random stream depending on the clock value (unpredictable).
myRand = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(myRand);

%% Dialog Box for User Input

prompt = {'Number of Permutations', 'Sampling Rate', '(m)ixed/(r)epeated', 'Number of Groups?' 'Results Name'};
dlg_title = 'Define Parameters';
def = {num2str(nPerm), num2str(rSample), type, num2str(nGroups), saveName};
noGr = inputdlg(prompt, dlg_title, 1, def);


nPerm      	= str2double(noGr{1});  % number of permutations 
rSample    	= str2double(noGr{2});  % sampling rate
type 		= noGr{3};
nGroups		= str2double(noGr{4});  % number of first factor levels (i.e groups) 
saveName   	= noGr{5};


% Data input
if nargin < 1;
    
    for i = 1:nGroups
        [DataFile{i}, DataPath{i}] = uigetfile('*.mat', ['Please Select Group ', num2str(i), ' Conditions'], 'MultiSelect', 'on');
    end
    
    if DataFile{1} == 0
        fprintf(1, 'Returning: No Data File Selected... \n');
        return;
    end
    
tic; % Start the timer for the entire analysis

    nLevels = cell2mat(cellfun(@(x) numel(x), DataFile, 'UniformOutput', false));
    
display ('Checking Input Files...')    
    if numel(unique(nLevels)) > 1
        error('All groups must have the same number of conditions');
    end
    
display ('Loading EEG Data...')
    try
        for j = 1:nGroups;
            for i = 1:length(DataFile{j})
                
                FullFileName    = strcat(DataPath{j}, DataFile{j}{i});
                LoadData        = load (FullFileName);
                Data{i,j}       = LoadData.Summary;
                DataLabel{i,j}  = DataFile{j}{i};
            end
        end
    catch
        error('Could not load the EEG files; possibly data format issue?')
    end
    
    [ElecFile, ElecPath] = uigetfile('', 'Please Select the Electrodes File');
    FullElecName = strcat(ElecPath, ElecFile);
    load (FullElecName);

else
    error('Script does not currently accept input arguments');
end

%% Rotate the data file so that groups are columns
Data = Data';
 
%% Compute Data Sizes

nP = cell2mat(cellfun(@(x) size(x,1), Data, 'UniformOutput', false)); % Find the number of participants in each dataset

display ('Checking Data...')
    if numel(unique(nP)) > 1
        error('All datsets must have the same number of participants (Balanced Design)');
    end

nP   = unique(nP);
nCh  = size(Data{1},2);
nS   = size(Data{1},3);


%% Error Checking
if ~isequal(nCh, length(e_loc))
    error ('Number of channels in data does not equal that of locations file')
end

%% Calculate the channels neighbours... using the modified version ChN2

display ('Calculating Channel Neighbours...')
ChN = ept_ChN2(e_loc);

%% Create all variables in loop at their maximum size to increase performance
maxTFCE.A  = zeros(nPerm,1);
maxTFCE.B  = zeros(nPerm,1);
maxTFCE.AB = zeros(nPerm,1);

%% Calculate the actual T values of all data
% Calculate different T-values for mixed and repeated measures ANOVA
display ('Calculating Observed Statistics...')
switch type
    case 'm'
       
        F_Obs = ept_mixedANOVA(Data);
        
    case 'r'

        F_Obs = ept_rmANOVA(Data);

    otherwise

        error('Analysis type must be (m)ixed or (r)epeated measures designs.\n')
end

% TFCE transformation...

TFCE_Obs.A  = ept_mex_TFCE2D(F_Obs.A,  ChN, E_H);
TFCE_Obs.B  = ept_mex_TFCE2D(F_Obs.B,  ChN, E_H);
TFCE_Obs.AB = ept_mex_TFCE2D(F_Obs.AB, ChN, E_H);

%% Calculating the T value and TFCE enhancement of each different permutation
display ('Calculating Permutations...')

    for i   = 1:nPerm
        
            switch type
                case 'm'
                    F_Perm = ept_mixedANOVA(Data, 1);
                case 'r'
                    F_Perm = ept_rmANOVA(Data, 1);
            end

        TFCE_Perm.A  = ept_mex_TFCE2D(F_Perm.A,  ChN, E_H);
        TFCE_Perm.B  = ept_mex_TFCE2D(F_Perm.B,  ChN, E_H);
        TFCE_Perm.AB = ept_mex_TFCE2D(F_Perm.AB, ChN, E_H);
        
        maxTFCE.A(i)  = max(max(TFCE_Perm.A));       % stores the maximum absolute value
        maxTFCE.B(i)  = max(max(TFCE_Perm.B)); 
        maxTFCE.AB(i) = max(max(TFCE_Perm.AB)); 
        
        progressbar(i/nPerm); %Progress Bar
       
    end

%% Calculating the p value from the permutation distribution
display ('Calculating Final Statistics...')

% add observed maximum
edges.A  = [maxTFCE.A;  max(max(TFCE_Obs.A))];
edges.B  = [maxTFCE.B;  max(max(TFCE_Obs.B))];
edges.AB = [maxTFCE.AB; max(max(TFCE_Obs.AB))];


[~,bin.A]      = histc(abs(TFCE_Obs.A),sort(edges.A));
P_Values.A     = 1-bin.A./(nPerm+2);
[~,bin.B]      = histc(abs(TFCE_Obs.B),sort(edges.B));
P_Values.B     = 1-bin.B./(nPerm+2);
[~,bin.AB]     = histc(abs(TFCE_Obs.AB),sort(edges.AB));
P_Values.AB    = 1-bin.AB./(nPerm+2);

% plot the "significance power"...
if plots == 1;
    figure; plot(sum(-log(P_Values.A)));  title ('Main Effect: Factor A')
    figure; plot(sum(-log(P_Values.B)));  title ('Main Effect: Factor B')
    figure; plot(sum(-log(P_Values.AB))); title ('Interaction Effect: AB')
end

% Save test information in single structure
display ('Saving Data...')
c = clock;
Info.Comments = ['TFCE analysis conducted at ', num2str(c(4)), ':', num2str(c(5)), ' on ', date];

Info.Parameters.E_H         = E_H;
Info.Parameters.nPerm       = nPerm;
Info.Parameters.rSample     = rSample;
Info.Parameters.type        = type;
Info.Parameters.nChannels   = nCh;
Info.Parameters.nSamples    = nS;

if size(Data,1)==2
    Info.Parameters.GroupSizes  = [size(Data{1,1}, 1); size(Data{2,1}, 1)];
elseif size(Data,1)==3
    Info.Parameters.GroupSizes  = [size(Data{1,1}, 1); size(Data{2,1}, 1); size(Data{3,1},1)];   
end

Info.Electrodes.e_loc               = e_loc;
Info.Electrodes.ChannelNeighbours   = ChN;

Info.DataFiles              = DataLabel;
Info.DataPaths              = DataPath;

Results.Obs                 = F_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = maxTFCE;
Results.P_Values            = P_Values;

save (saveName, 'Info', 'Data', 'Results', '-mat')

toc
%% Brief Report Summary in Matlab
%FactorA
[min_P, idx] = min(Results.P_Values.A(:));
[Ch, S]      = ind2sub(size(Results.P_Values.A),idx);
max_Obs      = Results.Obs.A(idx);

display(['FactorA peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);

% FactorB
[min_P, idx] = min(Results.P_Values.B(:));
[Ch, S]      = ind2sub(size(Results.P_Values.B),idx);
max_Obs      = Results.Obs.B(idx);

display(['FactorB peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,2)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);

% Interaction
[min_P, idx] = min(Results.P_Values.AB(:));
[Ch, S]      = ind2sub(size(Results.P_Values.AB),idx);
max_Obs      = Results.Obs.AB(idx);

display(['Interaction peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);


%% 

