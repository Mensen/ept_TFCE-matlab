function main_effects = ept_lme_main_effects (best_model, key_predictors)
% generates a table of main effects for all terms with key_predictors
% compares the model without the term to the "best model" using a Likelihood Ratio test

% valid only for two-way interactions

% get terms from best_model
model_term_names = best_model.Formula.FELinearFormula.TermNames;

% get full_table from the model
full_table = best_model.Variables;

% generate random_tag from model
random_factor_name = {best_model.Formula.GroupingVariableNames(:)};
random_factor_name = random_factor_name{1};

random_tag = ' ';
for n = 1 : length(random_factor_name)
    random_tag = [random_tag, ...
        ' + ( 1 | ' random_factor_name{n}{1}, ' )'];
end

% get r_squared of best_model
best_rsquared = 1 - best_model.SSE / best_model.SST;

% main effects of key predictors
% ''''''''''''''''''''''''''''''
main_effects = table();
for k = 1 : length(key_predictors)
    
    % find key predictor terms and remove them
    tmp_ind = contains(model_term_names, ...
        key_predictors{k});
    to_remove = model_term_names(tmp_ind);
    
    % remove terms
    new_fe_description = best_model.Formula.FELinearFormula;
    for n = 1 : length(to_remove)
        new_fe_description = removeTerms(new_fe_description, to_remove{n});
    end
    
    % calculate new model without key predictor{k}
    new_description = [new_fe_description.char, random_tag];
    new_model = fitlme(full_table, new_description);
    
    % calculate R-squared
    new_rsquared = 1 - new_model.SSE / new_model.SST;
    
    % compare to best_model
    LR_stats = compare(new_model, best_model);
    
    tmp_table = table();
    tmp_table.key_predictor = key_predictors(k);
    tmp_table.R2 = new_rsquared;
    tmp_table.delta_R2 = best_rsquared - new_rsquared;
    tmp_table.LRStat = LR_stats.LRStat(2);
    tmp_table.deltaDF = LR_stats.deltaDF(2);
    tmp_table.pValue = LR_stats.pValue(2);

    % add to group table
    main_effects = [main_effects; tmp_table];

end

% interaction effects of key predictors
% '''''''''''''''''''''''''''''''''''''
% list model terms with interactions between key_predictor

% find key predictor terms and remove them
tmp_ind = logical(contains(model_term_names, key_predictors) .* ...
    contains(model_term_names, ':'));
to_remove = model_term_names(tmp_ind);

for k = 1 : length(to_remove)

    % remove terms
    new_fe_description = best_model.Formula.FELinearFormula;
    new_fe_description = removeTerms(new_fe_description, to_remove{k});
    
        % calculate new model without key predictor{k}
    new_description = [new_fe_description.char, random_tag];
    new_model = fitlme(full_table, new_description);
    
    % calculate R-squared
    new_rsquared = 1 - new_model.SSE / new_model.SST;
    
    % compare to best_model
    LR_stats = compare(new_model, best_model);
    
    tmp_table = table();
    tmp_table.key_predictor = to_remove(k);
    tmp_table.R2 = new_rsquared;
    tmp_table.delta_R2 = best_rsquared - new_rsquared;    
    tmp_table.LRStat = LR_stats.LRStat(2);
    tmp_table.deltaDF = LR_stats.deltaDF(2);
    tmp_table.pValue = LR_stats.pValue(2);

    % append to group table
    main_effects = [main_effects; tmp_table];

end