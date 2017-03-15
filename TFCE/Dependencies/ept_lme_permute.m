function perm_table = ept_lme_permute(table, factors)
% permutes the input table according to factors indicated and type of 
% randomization permissible (between or within randomisation)

% TODO: current factors hard-coded for stroke experiment, need to generalize

% permutation of factor for all values for that participant
unique_participants = unique(table.participant);

% copy the table
perm_table = table;

% between factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
between_variables_names = {factors(find(~[factors.flag_within])).name};

% TODO: loop in case of multiple between subject factors

% find the unique group values
between_values = nan(length(unique_participants), 1);
for n = 1 : length(unique_participants)
   
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    % look at the group variable for that participant
    between_values(n) = unique(table.(between_variables_names{1})(rows));
    
end

% create random permutation of the group labels
% -1 because conversion from nominal adds 1 (*I think)
randomised_group = randsample(between_values - 1, length(between_values), 0);

% re-assign the labels to the participant in the table
for n = 1 : length(unique_participants)
    
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    % assign the new values to the between variable
    perm_table.(between_variables_names{1})(rows) = nominal(randomised_group(n));
      
end

% assign corresponding values to left/right groups
% TODO: figure out how to generalise the linking of between-factors
for n = 0 : 3
    indices = perm_table.(between_variables_names{1}) == num2str(n);

    if n == 1 || n == 3
        perm_table.left_stroke(indices) = '1';
    else
        perm_table.left_stroke(indices) = '0';
    end
    
    if n == 2 || n == 3
        perm_table.right_stroke(indices) = '1';
    else
        perm_table.right_stroke(indices) = '0';
    end

end
% eliminate ghost labels
perm_table.left_stroke = nominal(perm_table.left_stroke);
perm_table.right_stroke = nominal(perm_table.right_stroke);

% within factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
within_variables_names = {factors(find([factors.flag_within])).name};

for n = 1 : length(unique_participants)
    
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    randomised_trials = randsample(table.(within_variables_names{1})(rows), length(table.(within_variables_names{1})(rows)), 0);
    
    % assign the new value
    perm_table.(within_variables_names{1})(rows) = randomised_trials;   
    
end