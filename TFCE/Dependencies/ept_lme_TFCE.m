function p_values = ept_lme_TFCE(data, perm_data, channel_neighbours)
% run TFCE over observed data and each permutation

% get channel neighbourhood 
% observed_data
TFCE_Obs = ept_mex_TFCE2D(data, channel_neighbours, [0.66, 1]);

% permuted data
num_perm = size(perm_data, 1)
TFCE_max = nan(num_perm, 1);

for n = 1 : size(perm_data, 1)

    % enhance values using TFCE
    TFCE_perm = ept_mex_TFCE2D(perm_data(:, n), channel_neighbours, [0.66, 1]);

    % get the maximum value
    TFCE_max(n) = max(abs(TFCE_perm));
    
end

% calculate p-values
[~, bin] = histc(abs(TFCE_Obs), sort(TFCE_max));
p_values = 1 - bin ./ (num_perm);