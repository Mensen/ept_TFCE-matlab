%% EEG Permutation Test (egt) using Threshold-Free Cluster Enhancement 
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
% This tests initially computes a T-Value for each channel and sample
% These values are then enhanced or reduced depending on the size of
% the T-value and the size of the cluster the channel belongs to at 
% various thresholds...
%     
% TFCE is essentially the sum of the clusters to the power of 0.5 multiplied
% by the height of the current threshold squared
% 
% [Input]
% Participant data should be organised into two factors
%     e.g. Dataset 1 is controls while Dataset 2 is patients
% Electrode locations file created using eeglab is required to 
% % calculate channel neighbours
% 
% Analysis Types
% i = independent sample T-Test
% d = dependent (paired) sample T-Test
% o = one-sample T-test
% c = Pearson's correlation (r)
% 
% This Summary file should be a "Participants x Channels x Samples" variable
% - Participants should be in order of the group they belong to...
% - Channels must be in the same order as the corresponding electrodes file
% - Samples should be in chronological order
% - Values in the variable correspond to the mV averaged ERP data
% 
% [Output]
% - Info Structure
% - Data Structure
% - Results Structure
%

% Revision History
%
% Version 2.2
% 12.09.2013
% - Calls the new TFCE mex file which simultaneously looks for negative
% values (should be significantly faster then two separate calls)
%
% Version 2.1
% 11.09.2013
% - Included if statements for 2D and 3D TFCE differentiation
%
% 22.01.2013
% Corrected bug in correlation analysis (check for type 'r' instead of 'c')
% Eliminated the use of the shuffle function for correlation
% 
% 06.11.2012
% Included a one-sample test and some error handling
%
% 14.10.2012
% Includes Info/Data/Results Structure for better overview of test
% Includes the type 'c' for correlational analysis
% 23.12.2010
% - Performs two TFCE calculations, positive and negative values and then
%   puts them back together
%       - To prevent non existent clusters forming when taking abs values
%
% 03.01.2011
% - Calculated significance power and displays graph
%       - Significance power is the inverse of the p-value squared...
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
%
% 16.12.2011
% - Used -log of p-values to plot (more interpretable than power values)
% - Assigns defaults rather than prompting user for information
% 

 function []= ept_TFCE(DataFile1, DataFile2, ElecFile, varargin)

% assignin('base', 'varargin', varargin); 
%% Set Defaults
 
 E_H        = [0.66 2]; % default parameters of E and H
 nPerm      = 5000; % default number of permutations
 rSample    = 250; % default assumed sampling rate (not used in calculation but only for result viewer)
 saveName   = ['ept_Results_', date, '.mat'];
 type       = 'i'; % i = independent sample T-Test; d = dependent (paired) sample T-Test; c = Pearson's correlation
 plots      = 0; % change to '1' to show significance plots after calculation
  
% Set random stream depending on the clock value (unpredictable).
myRand = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(myRand);
 
 %% Process Arguments
 if nargin < 1;
    
    prompt = {'Number of Permutations', 'Sampling Rate', '(i)n(d)ependent / (c)orrelation?', 'Results Name', 'Plots? 0=No, 1=Yes'};
    dlg_title = 'Define Parameters';
    def = {num2str(nPerm), num2str(rSample), type, saveName, num2str(plots)};
    noGr = inputdlg(prompt, dlg_title,1,def);
    
    if isempty(noGr)
        error('Analysis cancelled')
    end
            
    nPerm      = str2double(noGr{1});  % number of permutations 
    rSample    = str2double(noGr{2});  % sampling rate
    type       = noGr{3};
    saveName   = noGr{4};
    plots      = str2double(noGr{5});
     
    if type == 'i' || type == 'd'
        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the first EEG Summary File', 'MultiSelect', 'off');
        [DataFile{2}, DataPath{2}] = uigetfile('', 'Please Select the second EEG Summary File', 'MultiSelect', 'off');
    elseif type == 'c'
        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the EEG Summary File', 'MultiSelect', 'off');
        [DataFile{2}, DataPath{2}] = uigetfile('', 'Please Select the Behavioural Summary File', 'MultiSelect', 'off');
    elseif type == 'o'
        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the EEG Summary File', 'MultiSelect', 'off');
    end
    
    if DataFile{1} == 0
        fprintf(1, 'Returning: No Data File Selected... \n');
        return;
    end
    
    FullFileName1 = strcat(DataPath{1}, DataFile{1});
    Data{1,1}     = load (FullFileName1);
    Data{1,1}     = Data{1}.Summary;
    
    % Error handling for non-onesample tests...
    if type ~= 'o'
        FullFileName2 = strcat(DataPath{2}, DataFile{2});
        Data{2,1}     = load (FullFileName2);
    else
        Data{2,1}     = zeros(size(Data{1}));
    end
    
    % Error handling for special behavioural data exception
    if type == 'i' || type == 'd'
        Data{2,1}	= Data{2}.Summary;
        aData = [Data{1};Data{2}];
    elseif type == 'c'
        Data{2}     = Data{2}.Behavioural;
    end
    
    % Load the electrode file...
    [ElecFile, ElecPath] = uigetfile('', 'Please Select the Electrodes File');
    FullElecName = strcat(ElecPath, ElecFile);
    load (FullElecName);
    
 elseif nargin > 2; 
    
    % Not a lot of error checking done here yet.
    DataFile{1} = DataFile1;
    DataFile{2} = DataFile2;
    
    Data{1}       = load (DataFile1, 'Summary');
    Data{1}       = Data{1}.Summary;
    Data{2}       = load (DataFile2);

    if type == 'i' || type == 'd'
        Data{2} 	= Data{2}.Summary;
        aData = [Data{1};Data{2}];
    elseif type == 'c'
        Data{2}     = Data{2}.Behavioural;
    end
    
    e_loc         = load (ElecFile);
    e_loc         = e_loc.e_loc;
 
 end

% Process Secondary Arguments
if nargin > 2
  if (round(nargin/2) == nargin/2)
    error('Even number of input arguments??')
  end
  for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~ischar(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    
    switch Param
        case 'e_h'
            E_H         = Value;
        case 'nperm'
            nPerm       = Value;
        case 'rsample'
            rSample     = Value;
        case 'savename'
            saveName    = Value;
        case 'type'
            type        = lower(Value);
        case 'plots'
            plots       = Value;
    end
  end
end

nA   = size(Data{1},1);
nB   = size(Data{2},1);
nCh  = size(Data{1},2);
nS   = size(Data{1},3);

% For Frequency-Time Data...
if ndims(Data{1})==4;
    nS = size(Data{1},4);
    nF = size(Data{1},3);
else
    nF = [];
end

%% -- Error Checking -- %%
if ~isequal(nB, nA) && type == 'r'
    error ('Must have the same number of participants as behavioural points...')
end

if ~isequal(nB, nA) && type == 'd'
    error ('Must have the same number of participants for paired comparisons...')
end

% Check Location File for Number of Channels
if ~isequal(nCh, length(e_loc))
    error ('Number of channels in data does not equal that of locations file')
end

%% Create Summary File

% Summary = [Data{1};  Data{2}];
tic; % Start the timer for the entire analysis

%% Calculate the channels neighbours... using the modified version ChN2

display('Calculating Channel Neighbours...')
ChN = ept_ChN2(e_loc);
display('Done')

%% Create all variables in loop at their maximum size to increase performance

maxTFCE = zeros(nPerm,1);

%% Calculate the actual T values of all data

display('Calculating Actual Differences...')

% Calculate different T-values for independent and dependent groups
    if type == 'i'
        
        T_Obs = (mean(Data{1})-mean(Data{2}))./sqrt(var(Data{1})/nA+var(Data{2})/nB);
        T_Obs = squeeze(T_Obs);   

    elseif type == 'd' || type == 'o'

        D    = Data{1}-Data{2};

        T_Obs = mean(D,1)./(std(D)./sqrt(nA));
        T_Obs = squeeze(T_Obs);
        
    elseif type == 'c'
        
        %Repmat B to match the Data sizes
        B2 = Data{2}(:, ones(nCh,1), ones(nS,1));

        n = size(Data{2},1);
        Exy = sum(B2.*Data{1});
        ExEy = sum(Data{2})*sum(Data{1});
        Ex_2 = sum(Data{2})^2;
        E_x2 = sum(Data{2}.^2);
        Ey_2 = sum(Data{1}.^2);
        E_y2 = sum(Data{1}.^2);

        T_Obs = squeeze((Exy-(ExEy./n)) * 1./ (sqrt(E_x2-(Ex_2/n))*sqrt(E_y2-(Ey_2./n))));
     
    end

% TFCE transformation...
if ismatrix(T_Obs);
    TFCE_Obs = ept_mex_TFCE2D(T_Obs, ChN, E_H);
end

if ndims(T_Obs) == 3;
    TFCE_Obs = ept_mex_TFCE3D(T_Obs, ChN, E_H);
end

display('Done')

%% Calculating the T value and TFCE enhancement of each different permutation

display('Calculating Permutations...')

    for i   = 1:nPerm           % Look into parfor for parallel computing
        
            if type == 'i' %two sample independent T-test
                r_perm      = randperm(nA+nB); % Consider using Shuffle mex here (50-85% faster)...
                
                if ismatrix(T_Obs);
                    nData       = aData(r_perm,:,:); 
                    sData{1}    = nData(1:nA,:,:); sData{2}= nData((nA+1):(nA+nB),:,:);
                else
                    nData       = aData(r_perm,:,:,:); 
                    sData{1}    = nData(1:nA,:,:,:); sData{2}= nData((nA+1):(nA+nB),:,:,:);
                end 
              
                T_Perm = (mean(sData{1})-mean(sData{2}))./sqrt(var(sData{1})/nA+var(sData{2})/nB);
                T_Perm = squeeze(T_Perm);   
                
            elseif type == 'd' || type == 'o' %one-sample T-test
                Signs =[-1,1];
                SignSwitch = randsample(Signs,size(D,1),'true')';
                if ismatrix(T_Obs);
                   SignSwitch = repmat(SignSwitch, [1 nCh nS]);
                else
                   SignSwitch = repmat(SignSwitch, [1 nCh nF nS]);
                end 
                        
                nData = SignSwitch.*D;
        
                T_Perm = mean(nData,1)./(std(nData)./sqrt(size(nData,1)));
                T_Perm = squeeze(T_Perm);
                
        
            elseif type == 'c' %correlation analysis
                %Repmat B to match the Data sizes
                Bp  = Data{2}(randperm(nB));
                B2p = Bp(:, ones(nCh,1), ones(nS,1));

                n    = size(Bp,1);
                Exy  = sum(B2p.*Data{1});
                ExEy = sum(Bp)*sum(Data{1});
                Ex_2 = sum(Bp)^2;
                E_x2 = sum(Bp.^2);
                Ey_2 = sum(Data{1}).^2;
                E_y2 = sum(Data{1}.^2);

                T_Perm = squeeze((Exy-(ExEy./n)) * 1./ (sqrt(E_x2-(Ex_2/n))*sqrt(E_y2-(Ey_2./n))));

            else
                error('Unrecognised analysis-type; see help file for valid inputs')
            end
            
        % TFCE transformation...
        if ismatrix(T_Perm);
            TFCE_Perm = ept_mex_TFCE2D(T_Perm, ChN, E_H);
        end

        if ndims(T_Perm) == 3;
            TFCE_Perm = ept_mex_TFCE3D(T_Perm, ChN, E_H);
        end
        
        maxTFCE(i) = max(abs(TFCE_Perm(:)));       % stores the maximum absolute value
        
        progressbar(i/nPerm); %Progress Bar
       
    end

display('Done')
    
%% Calculating the p value from the permutation distribution

display('Calculating P-Values and Saving...')

% add observed maximum
edges = [maxTFCE;max(abs(TFCE_Obs(:)))];

[~,bin]     = histc(abs(TFCE_Obs),sort(edges));
P_Values    = 1-bin./(nPerm+2);

% plot the "significance power"...

if plots == 1;
    if ismatrix(T_Perm);
        figure
        plot(sum(-log(P_Values)));
        title (SaveName)
    end
end

% Save test information in single structure
c = clock;
Info.Comments = ['TFCE analysis conducted at ', num2str(c(4)), ':', num2str(c(5)), ' on ', date];

Info.Parameters.E_H         = E_H;
Info.Parameters.nPerm       = nPerm;
Info.Parameters.rSample     = rSample;
Info.Parameters.type        = type;
Info.Parameters.nChannels   = nCh;
Info.Parameters.nSamples    = nS; % Not sure if actually used...
Info.Parameters.GroupSizes  = [nA, nB];

Info.Electrodes.e_loc       = e_loc;
Info.Electrodes.ChannelNeighbours = ChN;

Info.DataFiles              = DataFile;

Results.Obs                 = T_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = sort(maxTFCE);
Results.P_Values            = P_Values;

save (saveName, 'Info', 'Data', 'Results', '-mat')

%%
display('All done!')
toc

[min_P, idx] = min(Results.P_Values(:));
[Ch, S]      = ind2sub(size(Results.P_Values),idx);
max_Obs      = Results.Obs(idx);

display(['Peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);


end