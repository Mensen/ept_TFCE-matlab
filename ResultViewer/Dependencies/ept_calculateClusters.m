function [cluster_results] = ept_calculateClusters(Results, ChN, threshold)

%% -- Get necessary parameters -- %%
P_Values    = Results.P_Values;
TFCE_Obs    = Results.TFCE_Obs;
Obs_Values  = Results.Obs;
cluster_results = [];

%% -- Process Data -- %%
% For 2D datasets (Channel*Time || Channel*Frequency
if ndims(P_Values) == 2
    
    % Look for new clusters
    x = TFCE_Obs; % copy TFCE_Obs
    x(P_Values> threshold) = 0; % Threshold the data at alpha
    Cp = ept_ClusRes(x, ChN, 0.01); % Calculate Positive Clusters
    
    xn = x; % copy x
    xn(x>0)=0; % Only show negative x
    xn = abs(xn); % Make negative x's positive
    Cn = ept_ClusRes(xn, ChN, 0.01); % Calculate Negative Clusters
    
    C = Cn-Cp; % Combine the two
       
    b = unique(C); % How many different clusters are there?
    b(b==0)=[]; % Eliminate the 0 from being a unique cluster
    
    if numel(b)==0;
        display (['There are no clusters of significant data at the p = ' num2str(threshold) ' threshold']);
        return
    else
        
        cluster_results = struct();
        
        for n = 1:size(b,1)
            
            x = TFCE_Obs;
            x(C~=b(n))      = 0; %
            idPeak          = find(abs(x)==max(abs(x(:))));
            [PeakC, PeakS]  = ind2sub(size(x),idPeak);
            
            idSize          = find(C== b(n)); % find the rows and columns that are significant
            [SizeC, SizeS]  = ind2sub(size(C),idSize);
            
            cluster_results(n).channel_peak = PeakC(1);
            cluster_results(n).sample_peak = PeakS(1);
            cluster_results(n).max_t_value = Obs_Values(PeakC(1),PeakS(1));
            cluster_results(n).p_value_peak = P_Values(PeakC(1),PeakS(1));
            cluster_results(n).cluster_size = numel(idSize);
            cluster_results(n).unique_channels = numel(unique(SizeC));
            cluster_results(n).unique_samples = numel(unique(SizeS));
            cluster_results(n).sample_range = [num2str(min(SizeS)), ' - ', num2str(max(SizeS))];
            
        end
    end
    
else % for Time-Freqency Data...
    
    % Calculate clusters above p-value threshold
    data = TFCE_Obs; % copy TFCE_Obs
    data(P_Values>threshold) = 0; % Threshold the data at alpha
    C = ept_ClusRes3D(data, ChN, 0.01); % Calculate Negative Clusters    
    
    b = unique(C); % How many different clusters are there?
    b(b==0)=[]; % Eliminate the 0 from being a unique cluster
    
    if numel(b)==0;
        display (['There are no clusters of significant data at the p = ' num2str(threshold) ' threshold']);
        return
    end
    
    for n = 1:size(b,1);  
        x = TFCE_Obs;
        x(C~=b(n))      = 0; %
        idPeak          = find(abs(x)==max(abs(x(:))));
        [PeakC, PeakF, PeakS] = ind2sub(size(x),idPeak);
        
        idSize          = find(C== b(n)); % find the rows and columns that are significant
        [SizeC, SizeF, SizeS]=ind2sub(size(C),idSize);

        ClusRes{n,1}    = PeakC(1); % peak channel (just the first of many possible peak channels (but averaging may result in a channel in between two that is not significant)!
        ClusRes{n,2}    = PeakS(1);
        ClusRes{n,3}    = PeakF(1);
        ClusRes{n,4}    = Obs_Values(PeakC(1),PeakF(1),PeakS(1));
        ClusRes{n,5}    = P_Values(PeakC(1),PeakF(1),PeakS(1));                
        ClusRes{n,6}    = numel(idSize);
        ClusRes{n,7}    = numel(unique(SizeC));
        ClusRes{n,8}    = numel(unique(SizeS));
        ClusRes{n,9}    = numel(unique(SizeF));
        ClusRes{n,10}   = [num2str(min(SizeS)), ' - ', num2str(max(SizeS))];
        ClusRes{n,11}   = [num2str(min(SizeF)), ' - ', num2str(max(SizeF))];
        
    end    
end    