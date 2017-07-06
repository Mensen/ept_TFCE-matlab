% combine different permutation runs into a single result

% select the files
[filename, filepath] = uigetfile('.mat', 'multiSelect', 'on');

% load first file
load(filename{1});

% loop over files
for n = 2 : length(filename)
    
    % load each file
    loaded = load(filename{n});
    
    % check observed to ensure its same data and model
    if isequal(observed.beta, loaded.observed.beta)
        
        perm_values.t_value = [perm_values.t_value, loaded.perm_values.t_value];
        perm_values.tfce_value = [perm_values.tfce_value, loaded.perm_values.tfce_value];

        factors.num_perm = factors.num_perm + loaded.factors.num_perm;
        
    else
        
       error('observed beta values are not equal, probably different data');
        
    end
    
end

% recalculate observed p_values based on total permutations
% determine significance of each fixed factor from t-values
num_coeff = length(full_model.Coefficients);
num_chans = length(e_loc);

for n = 1 : num_coeff
    % calculate the inverse percentile (how many perms are larger than original)
    for nch = 1 : num_chans
        observed.max_pvalue(n, nch) = [sum(perm_values.t_value(n, :) ...
            > abs(observed.t_value(n, nch))) + 1] / [factors.num_perm + 1];
    end
end

% determine significance of each fixed factor from tfce-values
for n = 1 : num_coeff
    % calculate the inverse percentile (how many perms are larger than original)
    for nch = 1 : num_chans
        observed.tfce_pvalue(n, nch) = [sum(perm_values.tfce_value(n, :) ...
            > abs(observed.tfce_value(n, nch))) + 1] / [factors.num_perm + 1];
    end
end

save('ndeTFCE_10_11Hz_cov', ...
    'desired_frequency_bin', ...
    'e_loc', ...
    'factors', ...
    'full_model', ...
    'model_description', ...
    'observed', ...
    'perm_values');
