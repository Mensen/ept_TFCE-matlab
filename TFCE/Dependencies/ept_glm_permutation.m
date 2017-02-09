function [perm_data] = ept_glm_permutation(data, table, model_description)

num_perm = 500;
num_chan = size(data, 2);
num_coeff = 12;

perm_data = struct(...
    'ss_ratio', nan(num_chan, num_perm), ...
    'beta_values', nan(num_coeff, num_chan, num_perm), ...
    'se_values', nan(num_coeff, num_chan, num_perm));

swa_progress_indicator('initiate', 'number of permutations complete')
for n_perm = 1 : num_perm

    % update progress
    swa_progress_indicator('update', n_perm, num_perm);
    
    % get new table from old table
    new_table = ept_lme_permute(table);
    
    % run model once to get design matrix
    perm_model = fitlme(new_table, model_description);
    perm_design = designMatrix(perm_model, 'fixed');
    
    % run on all channels
    [perm_data.ss_ratio(:, n_perm), ...
        perm_data.beta_values(:, :, n_perm), ...
        perm_data.se_values(:, :, n_perm)] = ...
        ept_run_GLM(data, perm_design);
end