function [stats, p_value] = ept_quick_permutation(data1, data2, num_perm, measure)

% check inputs
if nargin < 4
    measure = 'difference';
end

if nargin < 3
    num_perm = 100;
end


% if k is low enough... then generate all possibilities
if length(data1) < 14
    possibilities = enumerate([1, -1], length(data1));
    num_perm = size(possibilities, 1);
else
    possibilities = [];
end


% calculate the observed difference
switch measure
    case 'difference'
        
        stats.observed = mean(data1) - mean(data2);
        
end

% loop for each permutation
for n = 1 : num_perm
    
    if ~isempty(possibilities)
        random_shuffle = possibilities(n, :)';
    else
        random_shuffle = datasample([1, -1], length(data1))';
    end
    
    % shuffle the data
    new_data = [data1 - data2] .* random_shuffle;

    
    % calculate the randomised difference
    switch measure
        case 'difference'
        
            stats.perm_values(n) = mean(new_data);

    end
end

% calculate the significance
stats.p_value = [sum(abs(stats.perm_values) ...
    >= abs(stats.observed)) + 1] / [num_perm + 1];

stats.num_perm = num_perm;
p_value = stats.p_value;


function possibilities = enumerate(x, n)

m = length(x);
X = cell(1, n);
[X{:}] = ndgrid(x);
X = X(end : -1 : 1);
y = cat(n+1, X{:});
possibilities = reshape(y, [m^n, n]);