function [ss_ratio, beta_values, se_values] = ept_run_GLM(data, design_matrix)
% function to compute the whole model F-value from a standard ept Data structure

% turn Data into list
if isa(data, 'cell')
    data = cell2mat(data(:));
end

% calculate beta values
beta_values = pinv(design_matrix) * data;

% calculate the standard error for each beta
sigma_squared = sum((data - design_matrix * beta_values).^2)./(size(design_matrix, 2)-size(design_matrix, 1));
chol_inv = inv(chol(design_matrix' * design_matrix));

% loop for each sigma
se_values = nan(size(beta_values));
for n = 1 : size(data, 2)
    covariance_matrix = abs(sigma_squared(n)) * (chol_inv*chol_inv');
    se_values(:, n) = sqrt(diag(covariance_matrix));
end

% predicted data with no error
data_hat = design_matrix * beta_values;

% sum of squares
ss_full = diag(data_hat' * data)' - sum(data, 1).^2 / size(data, 1);
ss_error = diag(data' * data - data_hat' * data)';

% ratio of explained vs remaining error
ss_ratio = ss_full ./ ss_error;


% get the residuals of the model
% residuals = data - data_hat;

% true f-value
% ~~~~~~~~~~~~
% true f-value depends on the degrees of freedom but permutation does not

% degrees of freedom
% dof_full = 2 - 1;
% dof_error = size(data, 1) - dof_full - 1;

% get the actual f_value
% f_value = [ss_full/dof_full] ./ [ss_error/dof_error];

% f_value = [ss_full/dof_full] ./ [ss_error/dof_error];

% R notation for finding individual beta significance
% # using direct calculations
% vY <- as.matrix(dfData[, -2])[, 5]                        # dependent variable
% mX <- cbind(constant = 1, as.matrix(dfData[, -2])[, -5])  # design matrix
% 
% vBeta <- solve(t(mX)%*%mX, t(mX)%*%vY)                    # coefficient estimates
% dSigmaSq <- sum((vY - mX%*%vBeta)^2)/(nrow(mX)-ncol(mX))  # estimate of sigma-squared
% mVarCovar <- dSigmaSq*chol2inv(chol(t(mX)%*%mX))          # variance covariance matrix
% vStdErr <- sqrt(diag(mVarCovar))                          # coeff. est. standard errors

