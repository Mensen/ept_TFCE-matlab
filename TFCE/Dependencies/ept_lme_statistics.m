function observed_statistics = ept_lme_permutation(full_table, model_description)

% TODO: missing output | patient indices | fixed stats

if nargin < 2
    model_description = ...
        'power ~ left_stroke * right_stroke * spindle_type + (left_stroke * right_stroke | participant)';
end

% ---------------------------------------------- %
% loop the model for each channel + permutations %
% ---------------------------------------------- %
% loop for each channel...
num_chans = size(output(1).topo_power, 1);
num_perms = 1999;

% pre-allocation (permutation statistics are n + 1)
% TODO: pre-define the size of the fixed_stats (number of coefficients)
observed_statistics = nan(size(fixed_stats, 1), num_chans, num_perms + 1);

% initiate progress meter
swa_progress_indicator('initiate', 'number of channels complete')

for nCh = 1 : num_chans
    % update progress
    swa_progress_indicator('update', nCh, num_chans);

    % get the data for that channel
    data_output = {output(~patient_indices | pat_t1_indices).topo_power};
    data_output = cellfun(@(x) x(nCh, :), data_output, 'uniform', false);
    
    % use log10 to make data normally distributed
    dependent_variable = log10([cell2mat(data_output)']);
    
    full_table.power = dependent_variable;
    
    % run the actual model
    observed_model = fitlme(full_table, model_description);
    [~, ~, fixed_stats] = fixedEffects(observed_model);
    
    % extract T statistics, permute, and rerun model
    observed_statistics(:, nCh, 1) = fixed_stats(:, 4);
    
    for n = 1 : num_perms        
        % get the permuted table
        % TODO generalise the table permutation function
        new_table = ept_lme_permute(full_table, [], []);
        
        % run the new model
        perm_model = fitlme(new_table, model_description);
        
        % get the fixed effects
        [~, ~, perm_effects] = fixedEffects(perm_model);
        
        % save the T-values from each permutation
        observed_statistics(:, nCh, n + 1) = perm_effects(:, 4);
        
    end
end

% hist(perm_stats(10, :));