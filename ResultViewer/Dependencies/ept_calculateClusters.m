function [cluster_results] = ept_calculateClusters(Results, ChN, threshold)

% -- Get necessary parameters -- %%
P_Values    = Results.P_Values;
TFCE_Obs    = Results.TFCE_Obs;
Obs_Values  = Results.Obs;
cluster_results = [];

% -- Process Data -- %%
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
    
    all_clusters = Cn-Cp; % Combine the two
       
    unique_clusters = unique(all_clusters); % How many different clusters are there?
    unique_clusters(unique_clusters==0) = []; % Eliminate the 0 from being a unique cluster
    
    if numel(unique_clusters) == 0;
        fprintf('There are no clusters of significant data at the p = %0.3f threshold', threshold);
        return
        
    else
        % pre-allocate as struct
        cluster_results = struct();
        
        % loop for each cluster
        for n = 1 : size(unique_clusters,1)
            
            temp_observed = TFCE_Obs;
            temp_observed(all_clusters ~= unique_clusters(n)) = 0;
            % find the largest TFCE value
            [~, peak_ind] = max(abs(temp_observed(:)));
            [peak_channel, peak_sample] = ...
                ind2sub(size(temp_observed), peak_ind);
            
            % find the rows and columns that are significant
            sig_ind = find(all_clusters == unique_clusters(n)); 
            [sig_channels, sig_samples] = ind2sub(size(all_clusters),sig_ind);
            
            % assign to structure
            % TODO: currently if tied, choose first... could report all
            cluster_results(n).channel_peak = peak_channel(1);
            cluster_results(n).sample_peak = peak_sample(1);
            cluster_results(n).max_t_value = ...
                Obs_Values(peak_channel(1), peak_sample(1));
            cluster_results(n).p_value_peak = ...
                P_Values(peak_channel(1), peak_sample(1));
            cluster_results(n).cluster_size = numel(sig_ind);
            cluster_results(n).unique_channels = numel(unique(sig_channels));
            cluster_results(n).unique_samples = numel(unique(sig_samples));
            cluster_results(n).sample_range = ...
                [num2str(min(sig_samples)), ' - ', num2str(max(sig_samples))];
            
            % export the actual supra-cluster points
            cluster_results(n).cluster_locations = false(size(Results.P_Values));
            cluster_results(n).cluster_locations(...
                all_clusters == unique_clusters(n)) = true;
            
        end
    end
    
else % for Time-Freqency Data...
    
    % Calculate clusters above p-value threshold
    data = TFCE_Obs; % copy TFCE_Obs
    data(P_Values>threshold) = 0; % Threshold the data at alpha
    all_clusters = ept_ClusRes3D(data, ChN, 0.01); % Calculate Negative Clusters    
    
    unique_clusters = unique(all_clusters); % How many different clusters are there?
    unique_clusters(unique_clusters==0)=[]; % Eliminate the 0 from being a unique cluster
    
    if numel(unique_clusters)==0;
        display (['There are no clusters of significant data at the p = ' num2str(threshold) ' threshold']);
        return
    end
    
    for n = 1:size(unique_clusters,1);  
        x = TFCE_Obs;
        x(all_clusters~=unique_clusters(n))      = 0; %
        idPeak          = find(abs(x)==max(abs(x(:))));
        [peak_channel, PeakF, peak_sample] = ind2sub(size(x),idPeak);
        
        sig_ind          = find(all_clusters== unique_clusters(n)); % find the rows and columns that are significant
        [SizeC, SizeF, sig_samples]=ind2sub(size(all_clusters),sig_ind);

        ClusRes{n,1}    = peak_channel(1); % peak channel (just the first of many possible peak channels (but averaging may result in a channel in between two that is not significant)!
        ClusRes{n,2}    = peak_sample(1);
        ClusRes{n,3}    = PeakF(1);
        ClusRes{n,4}    = Obs_Values(peak_channel(1),PeakF(1),peak_sample(1));
        ClusRes{n,5}    = P_Values(peak_channel(1),PeakF(1),peak_sample(1));                
        ClusRes{n,6}    = numel(sig_ind);
        ClusRes{n,7}    = numel(unique(SizeC));
        ClusRes{n,8}    = numel(unique(sig_samples));
        ClusRes{n,9}    = numel(unique(SizeF));
        ClusRes{n,10}   = [num2str(min(sig_samples)), ' - ', num2str(max(sig_samples))];
        ClusRes{n,11}   = [num2str(min(SizeF)), ' - ', num2str(max(SizeF))];
        
    end    
end    