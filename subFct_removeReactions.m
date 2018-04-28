function [N, v0, r_names_left, param_intervals] = subFct_removeReactions(r_exclude, r_include, N, v0, param_intervals, r_names, verbose)
% This function reduces the model by removing a specified number of reactions.
% Reactions can be specified by the input arguments 'r_exclude' (all reactions in 'r_exclude' will be removed) or
% in 'r_include' (all reactions NOT in 'r_include' will be removed).
% Note that only one of these arguments can be specified at a time, the other argument must be an empty cell or matrix.
% Specifying both argument 'r_exclude' and 'r_include' may result in errors.
% Possible formats for 'r_exclude' and 'r_include' are:
%   1. Cellstring with reaction names
%   2. Boolean array indicating the positions to be excluded/included by 1s
%   3. Numerical array containing the positions as integer values.
% All other input arguments have the same format as in the main function skm().

%% Determine the format in which included/excluded reactions are given
r_include_cell_names        = 0;
r_include_boolean_indices   = 0;
r_include_numeric_indices   = 0;
r_exclude_cell_names        = 0;
r_exclude_boolean_indices   = 0;
r_exclude_numeric_indices   = 0;

if ~isempty(r_include)
    % Check if input argument r_include is a cell of strings, a boolean array or an array of indices:
    if iscellstr(r_include)
        r_include_cell_names            = 1;
        % Check whether r_include is a matrix:
    elseif length(size(r_include))<=2 && all(size(r_include)>=1)
        if length(r_include)==length(r_names) && min(r_include)>=0 && max(r_include)<=1
            r_include_boolean_indices   = 1;
        elseif max(r_include) <= length(r_names)
            r_include_numeric_indices   = 1;
        end
    end
elseif ~isempty(r_exclude)
    % Check if input argument r_exclude is a cell of strings, a boolean array or an array of indices:
    if iscellstr(r_exclude)
        r_exclude_cell_names            = 1;
        % Check whether r_exclude is a matrix:
    elseif length(size(r_exclude))<=2 && all(size(r_exclude)>=1)
        if length(r_exclude)==length(r_names) && min(r_exclude)>=0 && max(r_exclude)<=1
            r_exclude_boolean_indices   = 1;
        elseif max(r_exclude) <= length(r_names)
            r_exclude_numeric_indices   = 1;
        end
    end
else
    % If neither r_include nor r_exclude are given, don't remove any reactions:
    r_remove_index = false(1, length(r_names));
end


%% Determine which reactions should be removed and store indices in a boolean array
if r_include_cell_names
    r_remove_index = ~ismember(r_names, r_include);
    unknown_names = ~ismember(r_include, r_names);
    if verbose
        if sum(unknown_names)
            fprintf('The following reactions are not known:\n');
            disp(char(r_include(unknown_names)));
            fprintf('\n')
        end
    end
elseif r_include_boolean_indices
    r_remove_index = ~r_include;
elseif r_include_numeric_indices
    r_remove_index = true(size(r_names));
    r_remove_index(r_include) = 0;
elseif r_exclude_cell_names
    r_remove_index = ismember(r_names, r_exclude);
    unknown_names = ~ismember(r_exclude, r_names);
    if verbose
        if sum(unknown_names)
            fprintf('The following reactions are not known:\n');
            disp(char(r_exclude(unknown_names)));
            fprintf('\n')
        end
    end
elseif r_exclude_boolean_indices
    r_remove_index = r_exclude;
elseif r_exclude_numeric_indices
    r_remove_index = false(size(r_names));
    r_remove_index(r_exclude) = 1;
end

if sum(r_remove_index) == length(r_names)
    fprintf('User required that all reactions should be removed. The SK-model would then be empty!\nAll reactions will be left in the model instead.\n\n')
    r_remove_index = false(size(r_names));
end

if verbose
    if sum(r_remove_index)
        fprintf('Removing the following reactions:\n');
        disp(char(r_names(r_remove_index)));
        fprintf('\nWith index:\n');
        disp(num2str(find(r_remove_index)));
        fprintf('\n');
    end
end

%% Remove reactions from the model
N(:, r_remove_index)    = [];
v0(r_remove_index)      = [];
r_names_left = r_names(~r_remove_index);

param_types = fieldnames(param_intervals);
for i = 1:length(param_types)
    current_param_matrix = param_intervals.(param_types{i});
    current_param_matrix(r_remove_index,:) = [];
    param_intervals.(param_types{i}) = current_param_matrix;
end

