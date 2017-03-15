function [observed, perm_values] = ept_slow_lme_permutation(data_of_interest, full_table, model_description, factors, channel_neighbours)
% script to run the complete lme analysis for all channels etc

NUM_PERMS = 1000;
FLAG_TFCE = true;

% run full model on all channels
% ''''''''''''''''''''''''''''''
% run single full_model for parameters
temp_model = fitlme(full_table, model_description);

% pre-allocate model (TODO: find number of coefficients before)
num_coeff = length(temp_model.CoefficientNames);
num_chans = size(data_of_interest, 1);
num_observations = size(data_of_interest, 2);

% get dependent variable name from model description
dv_name = strtok(model_description, ' ');

% pre-allocate
observed = struct(...
    'beta', nan(num_coeff, num_chans), ...
    'se', nan(num_coeff, num_chans), ...
    't_value', nan(num_coeff, num_chans), ...
    'tfce_value', nan(num_coeff, num_chans), ...
    'model_pvalue', nan(num_coeff, num_chans), ...
    'max_pvalue', nan(num_coeff, num_chans), ...
    'tfce_pvalue', nan(num_coeff, num_chans));
    
% loop the model for each channel
fprintf(1, '\nTesting observed model...');
swa_progress_indicator('initiate', 'number of channels complete')
for nCh = 1 : num_chans
    swa_progress_indicator('update', nCh, num_chans);
    
    % put the channel in the table
    full_table.(dv_name) = data_of_interest(nCh, :)';
    
    % run the actual model
    full_model_object = fitlme(full_table, model_description);
    
    % extract parameters of interest
    observed.beta(:, nCh) = double(full_model_object.Coefficients(:, 2));
    observed.se(:, nCh) = double(full_model_object.Coefficients(:, 3));
    observed.t_value(:, nCh) = double(full_model_object.Coefficients(:, 4));
    observed.model_pvalue(:, nCh) = double(full_model_object.Coefficients(:, 6));
end

% run TFCE analysis over observed t-values
if FLAG_TFCE
    fprintf(1, '\nRunning TFCE analysis...\n');
    for n = 1 : num_coeff
        temp_value = ept_mex_TFCE2D(repmat(observed.t_value(n, :), [2, 1])', channel_neighbours, [0.66, 1]);
        observed.tfce_value(n, :) = temp_value(:, 1);
    end
end

% run permutations
% ''''''''''''''''
if NUM_PERMS > 1
   
    % pre-allocate empirical distribution
    perm_values = struct(...
        't_value', nan(num_coeff, NUM_PERMS), ...
        'tfce_value', nan(num_coeff, NUM_PERMS));
    perm_tvalue = nan(num_coeff, num_chans);
    perm_tfce = nan(num_coeff, num_chans);
    
    % initiate progress meter
    swa_progress_indicator('initiate', 'number of permutations complete')
    
    for n_perm = 1 : NUM_PERMS
        % update progress
        swa_progress_indicator('update', n_perm, NUM_PERMS);
        
        % permute the table randomly (but appropriately)
        perm_table = ept_lme_permute(full_table, factors);
        
        for nCh = 1 : num_chans
            % put the channel in the table
            perm_table.(dv_name) = data_of_interest(nCh, :)';
            
            % run the actual model
            full_model_object = fitlme(perm_table, model_description);
            
            % extract parameters of interest
            perm_tvalue(:, nCh) = double(full_model_object.Coefficients(:, 4));
        end
        
        % run TFCE analysis over permuted data
        if FLAG_TFCE
            for n = 1 : num_coeff
                temp_value = ept_mex_TFCE2D(repmat(perm_tvalue(n, :), [2, 1])', channel_neighbours, [0.66, 1]);
                perm_tfce(n, :) = temp_value(:, 1);
            end
        end
        
        % find the maximum t-value (for each parameter separately)
        perm_values.t_value(:, n_perm) = max(perm_tvalue, [], 2);
        perm_values.tfce_value(:, n_perm) = max(perm_tfce, [], 2);
        
    end
end

% determine significance of each fixed factor from t-values
for n = 1 : num_coeff
    [~, bin] = histc(abs(observed.t_value(n, :)), sort(abs(perm_values.t_value(n, :))));
    observed.max_pvalue(n, :) = [bin + 1] ./ (NUM_PERMS + 1);
end

% determine significance of each fixed factor from t-values
if FLAG_TFCE
    for n = 1 : num_coeff
        [~, bin] = histc(abs(observed.tfce_value(n, :)), sort(abs(perm_values.tfce_value(n, :))));
        observed.tfce_pvalue(n, :) = [bin + 1] ./ (NUM_PERMS);
    end
end