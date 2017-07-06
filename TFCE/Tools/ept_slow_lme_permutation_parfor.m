function [observed, perm_values] = ept_slow_lme_permutation(data_of_interest, full_table, model_description, factors, channel_neighbours)
% script to run the complete lme analysis for all channels etc

% define defaults
NUM_PERMS = 0;
FLAG_TFCE = true;
FLAG_RANDOM = true;

% check arguments
if nargin < 5
    FLAG_TFCE = false; 
    fprintf(1, 'Warning: without the channels neighbourhood you cannot run TFCE');
end

% check for NUM_PERMS user input
if nargin > 3 & isfield(factors, 'num_perm')
    NUM_PERMS = factors.num_perm;
end

% assign output to perm_values in case none are run 
perm_values = [];

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

observed_beta = nan(num_coeff, num_chans);
observed_se = nan(num_coeff, num_chans);
observed_t_value = nan(num_coeff, num_chans);
observed_model_pvalue = nan(num_coeff, num_chans);

% loop the model for each channel
fprintf(1, '\nTesting observed model...');
parfor nCh = 1 : num_chans
    
    fprintf(1, '.\n');
    
    % put the channel in the table (required for parfor not to call struct directly
    temp_table = full_table;
    temp_table.(dv_name) = data_of_interest(nCh, :)';
    
    % run the actual model
    if FLAG_RANDOM
        full_model_object = fitlme(temp_table, model_description, ...
            'covariancePattern', 'Diagonal'); % FullCholesky | CompSymm | Diagonal | Full        
    else
        full_model_object = fitlme(temp_table, model_description);
    end
    
    % extract parameters of interest
    observed_beta(:, nCh) = double(full_model_object.Coefficients(:, 2));
    observed_se(:, nCh) = double(full_model_object.Coefficients(:, 3));
    observed_t_value(:, nCh) = double(full_model_object.Coefficients(:, 4));
    observed_model_pvalue(:, nCh) = double(full_model_object.Coefficients(:, 6));
end
fprintf(1, ' complete\n');

% for parfor put the observed values into the observed array
observed = struct(...
    'beta', observed_beta, ...
    'se', observed_se, ...
    't_value', observed_t_value, ...
    'model_pvalue', observed_model_pvalue);

% run TFCE analysis over observed t-values
if FLAG_TFCE
    fprintf(1, '\nRunning TFCE analysis...\n');
    for n = 1 : num_coeff
        temp_value = ept_mex_TFCE2D(repmat(observed.t_value(n, :), [2, 1])', channel_neighbours, [0.66, 1]);
        observed.tfce_value(n, :) = temp_value(:, 1);
    end
end

% pre-allocate fail detector
failed_perm = false(NUM_PERMS, 1);

% run permutations
% ''''''''''''''''
if NUM_PERMS > 0
   
   
    % pre-allocate empirical distribution
    perm_values = struct(...
        't_value', nan(num_coeff, NUM_PERMS), ...
        'tfce_value', nan(num_coeff, NUM_PERMS));
    perm_tvalue = nan(num_coeff, num_chans);
    perm_tfce = nan(num_coeff, num_chans);

    % start basic save file for temporary saving
    save('ept_temp_perm_save.mat', 'observed', 'perm_values', 'model_description', 'factors');
    
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
            
            % some permutations seem to be weird so...
            try
                % run the actual model
                full_model_object = fitlme(perm_table, model_description, ...
                    'covariancePattern', 'Diagonal'); % FullCholesky | CompSymm | Diagonal | Full
                
                % extract parameters of interest
                perm_tvalue(:, nCh) = double(full_model_object.Coefficients(:, 4));
            
            catch
                % if the model fails...
                failed_perm(n_perm) = true;
                % skip the rest of the channels since TFCE won't work anyway
                break 
            end
        end
        
        if ~failed_perm(n_perm) % check for failed permutation
            % run TFCE analysis over permuted data
            if FLAG_TFCE
                for n = 1 : num_coeff
                    temp_value = ept_mex_TFCE2D(repmat(perm_tvalue(n, :), [2, 1])', channel_neighbours, [0.66, 1]);
                    perm_tfce(n, :) = temp_value(:, 1);
                end
            end
            
            % find the maximum t-value (for each parameter separately)
            perm_values.t_value(:, n_perm) = max(abs(perm_tvalue), [], 2);
            perm_values.tfce_value(:, n_perm) = max(abs(perm_tfce), [], 2);
        end
        
        % update temp save with the permutation data in case of loss
        save('ept_temp_perm_save.mat', 'perm_values', '-append');
        
    end
else
    return
end

% TODO: eliminate the nan runs
nan_runs = isnan(perm_values.t_value(1, :));
perm_values.t_value(:, nan_runs) = [];
perm_values.tfce_value(:, nan_runs) = [];

% determine significance of each fixed factor from t-values
for n = 1 : num_coeff
    % calculate the inverse percentile (how many perms are larger than original)
    for nch = 1 : num_chans
        observed.max_pvalue(n, nch) = [sum(perm_values.t_value(n, :) ...
            > abs(observed.t_value(n, nch))) + 1] / [NUM_PERMS + 1];
    end
end

% determine significance of each fixed factor from tfce-values
if FLAG_TFCE
    for n = 1 : num_coeff
        % calculate the inverse percentile (how many perms are larger than original)
        for nch = 1 : num_chans
            observed.tfce_pvalue(n, nch) = [sum(perm_values.tfce_value(n, :) ...
                > abs(observed.tfce_value(n, nch))) + 1] / [NUM_PERMS + 1];
        end
    end
end


%         [~, ~, bin] = histc(abs(observed.tfce_value(n, :)), ...
%             sort(abs(perm_values.tfce_value(n, :))));
%         observed.tfce_pvalue(n, :) = 1 -[bin + 1] ./ (NUM_PERMS);
