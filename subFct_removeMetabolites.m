function [N, S0, m_dependent, m_names_left, param_intervals] = subFct_removeMetabolites(m_exclude, m_include, N, S0, m_dependent, param_intervals, m_names, verbose)
% This function reduces the model by removing a specified number of metabolites.
% Metabolites can be specified by the input arguments 'm_exclude' (all metabolites in 'm_exclude' will be removed) or
% in 'm_include' (all metabolites NOT in 'm_include' will be removed).
% Note that only one of these arguments can be specified at a time, the other argument must be an empty cell or matrix.
% Specifying both argument 'm_exclude' and 'm_include' may result in errors.
% Possible formats for 'm_exclude' and 'm_include' are:
%   1. Cellstring with metabolite names
%   2. Boolean array indicating the positions to be excluded/included by 1s
%   3. Numerical array containing the positions as integer values.
% All other input arguments have the same format as in the main function skm().


%% Determine the format in which included/excluded metabolites are given
m_include_cell_names        = 0;
m_include_boolean_indices   = 0;
m_include_numeric_indices   = 0;
m_exclude_cell_names        = 0;
m_exclude_boolean_indices   = 0;
m_exclude_numeric_indices   = 0;

if ~isempty(m_include)
    % Check if input argument m_include is a cell of strings, a boolean array or an array of indices:
    if iscellstr(m_include)
        m_include_cell_names            = 1;
        % Check whether m_include is a matrix:
    elseif length(size(m_include))<=2 && all(size(m_include)>=1)
        if length(m_include)==length(m_names) && min(m_include)>=0 && max(m_include)<=1
            m_include_boolean_indices   = 1;
        elseif max(m_include) <= length(m_names)
            m_include_numeric_indices   = 1;
        end
    end
elseif ~isempty(m_exclude)
    % Check if input argument m_exclude is a cell of strings, a boolean array or an array of indices:
    if iscellstr(m_exclude)
        m_exclude_cell_names            = 1;
        % Check whether m_exclude is a matrix:
    elseif length(size(m_exclude))<=2 && all(size(m_exclude)>=1)
        if length(m_exclude)==length(m_names) && min(m_exclude)>=0 && max(m_exclude)<=1
            m_exclude_boolean_indices   = 1;
        elseif max(m_exclude) <= length(m_names)
            m_exclude_numeric_indices   = 1;
        end
    end
else
    % If neither m_include nor m_exclude are given, don't remove any metabolites:
    m_remove_index = false(1, length(m_names));
end

%% Determine which metabolites should be removed and store indices in a boolean array
if m_include_cell_names
    m_remove_index = ~ismember(m_names, m_include);
    unknown_names = ~ismember(m_include, m_names);
    if verbose
        if sum(unknown_names)
            fprintf('The following metabolites are not known:\n');
            disp(char(m_include(unknown_names)));
            fprintf('\n')
        end
    end
elseif m_include_boolean_indices
    m_remove_index = ~m_include;
elseif m_include_numeric_indices
    m_remove_index = true(size(m_names));
    m_remove_index(m_include) = 0;
elseif m_exclude_cell_names
    m_remove_index = ismember(m_names, m_exclude);
    unknown_names = ~ismember(m_exclude, m_names);
    if verbose
        if sum(unknown_names)
            fprintf('The following metabolites are not known:\n');
            disp(char(m_exclude(unknown_names)));
            fprintf('\n')
        end
    end
elseif m_exclude_boolean_indices
    m_remove_index = m_exclude;
elseif m_exclude_numeric_indices
    m_remove_index = false(size(m_names));
    m_remove_index(m_exclude) = 1;
end

if sum(m_remove_index) == length(m_names)
    fprintf('User required that all metabolites should be removed. The SK-model would then be empty!\nAll metabolites will be left in the model instead.\n\n')
    m_remove_index = false(size(m_names));
end

if verbose
    if sum(m_remove_index)
        fprintf('Removing the following metabolites:\n');
        disp(char(m_names(m_remove_index)));
        fprintf('\nWith index:\n');
        disp(num2str(find(m_remove_index)));
        fprintf('\n');
    end
end

%% Remove metabolites from the model
N(m_remove_index, :)        = [];
S0(m_remove_index)          = [];
m_dependent(m_remove_index) = [];
m_names_left = m_names(~m_remove_index);

param_types = fieldnames(param_intervals);
for i = 1:length(param_types)
    current_param_matrix = param_intervals.(param_types{i});
    current_param_matrix(:,m_remove_index) = [];
    param_intervals.(param_types{i}) = current_param_matrix;
end