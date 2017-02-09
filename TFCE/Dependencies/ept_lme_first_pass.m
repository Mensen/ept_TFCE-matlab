function [full_model, reduced_data, design_fixed, full_model_object] = ...
    ept_lme_first_pass(data_of_interest, full_table, model_description)
% function to run the full model and reduced model once per input (channel) and
% return the model parameters for both

% check inputs
if nargin < 3
    model_description{1} = ...
        'power ~ left_stroke * right_stroke * spindle_type + (left_stroke * right_stroke | participant)';
    model_description{2} = ...
        'power ~ left_stroke * right_stroke * spindle_type';
end

% pre-allocate model (TODO: find number of coefficients before)
num_coeff = 12;
num_chans = size(data_of_interest{1}, 1);
num_observations = sum(cellfun(@(x) size(x, 2), data_of_interest));

full_model = struct(...
    'beta', nan(num_coeff, num_chans), ...
    'se', nan(num_coeff, num_chans));

reduced_data = nan(num_observations, num_chans);

% loop the model for each channel
% initiate progress meter
swa_progress_indicator('initiate', 'number of channels complete')

for nCh = 1 : num_chans
    % update progress
    swa_progress_indicator('update', nCh, num_chans);

    % get the data for that channel
    data_output = cellfun(@(x) x(nCh, :), data_of_interest, 'uniform', false);
    
    % use log10 to make data normally distributed
    dependent_variable = log10([cell2mat(data_output)']);
    full_table.power = dependent_variable;
    
    % run the actual model
    full_model_object = fitlme(full_table, model_description);
    
    % extract parameters of interest
    full_model.beta(:, nCh) = double(full_model_object.Coefficients(:, 2));
    full_model.se(:, nCh) = double(full_model_object.Coefficients(:, 3));
    
    % run the reduced model
    design_random = designMatrix(full_model_object, 'random');
    
    % get the residuals after random effects have been removed
    [random_beta, ~, ~] = randomEffects(full_model_object);
    nuisance_predicted = design_random * random_beta;
    reduced_data(:, nCh) = full_table.power - nuisance_predicted;

end

% get the design matrix for fixed and random effects
design_fixed = designMatrix(full_model_object, 'fixed');

