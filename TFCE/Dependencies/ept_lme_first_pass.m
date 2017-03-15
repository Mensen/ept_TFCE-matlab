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

% run single full_model for parameters
temp_model = fitlme(full_table, model_description);

% pre-allocate model (TODO: find number of coefficients before)
num_coeff = length(temp_model.CoefficientNames);
num_chans = size(data_of_interest, 1);
num_observations = size(data_of_interest, 2);

full_model = struct(...
    'beta', nan(num_coeff, num_chans), ...
    'se', nan(num_coeff, num_chans));

reduced_data = nan(num_observations, num_chans);

% get dependent variable name from model description
dv_name = strtok(model_description, ' ');

% loop the model for each channel
% initiate progress meter
swa_progress_indicator('initiate', 'number of channels complete')

for nCh = 1 : num_chans
    % update progress
    swa_progress_indicator('update', nCh, num_chans);
   
    % put the channel in the table
    full_table.(dv_name) = data_of_interest(nCh, :)';
    
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

