function new_table = ept_lme_permute(table, factors, type)
% permutes the input table according to factors indicated and type of 
% randomization permissible (between or within randomisation)

% TODO: current factors hard-coded for stroke experiment, need to generalize

% permutation of factor for all values for that participant
unique_participants = unique(table.participant);

% between factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% find the unique group values
group_values = nan(length(unique_participants), 1);
for n = 1 : length(unique_participants)
   
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    % look at the group variable for that participant
    group_values(n) = unique(table.stroke_location(rows));
    
end

% create random permutation of the group labels
% -1 because conversion from nominal adds 1 (*I think)
randomised_group = randsample(group_values - 1, length(group_values), 0);

% re-assign the labels to the participant in the table
new_table = table;
for n = 1 : length(unique_participants)
    
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    % assign the new values to stroke location
    new_table.stroke_location(rows) = nominal(randomised_group(n));
      
end

% assign corresponding values to left/right groups
for n = 0 : 3
    indices = new_table.stroke_location == num2str(n);

    if n == 1 || n == 3
        new_table.left_stroke(indices) = true;
    else
        new_table.left_stroke(indices) = false;
    end
    
    if n == 2 || n == 3
        new_table.right_stroke(indices) = true;
    else
        new_table.right_stroke(indices) = false;
    end

    
end


% within factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
for n = 1 : length(unique_participants)
    
    % get the relevant rows
    rows = table.participant==unique_participants(n);
    
    randomised_trials = randsample(table.spindle_type(rows), length(table.spindle_type(rows)), 0);
    
    % assign the new value
    new_table.spindle_type(rows) = randomised_trials;   
    
end