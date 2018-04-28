function [values_sorted, names_sorted] = subFct_concatenateParams(modelParam_values, modelParam_names)
% Concatenate feature names to create cell arrays.
% This function is invoked when generating classifier input and for performing pairwise tests on the parameters.

% Concatenate all parameter names:
names_reduced = modelParam_names;
param_types = fieldnames(names_reduced);
for type_number = 1:length(param_types)
    current_type    = param_types{type_number};
    current_field   = names_reduced.(current_type);
    if isempty(current_field)
        names_reduced = rmfield(names_reduced, current_type);
    end
end

names_cell = {};
param_types = fieldnames(names_reduced);
for type_number = 1:length(param_types)
    current_type    = param_types{type_number};
    names_cell = [names_cell; names_reduced.(current_type)];
end

if ~isempty(names_cell)
    names_str = strcat(names_cell(:,1), '-',names_cell(:,2), '-',names_cell(:,3));
else
    names_str = cell(0);
end

% Concatenate all values:
values_array = [modelParam_values.Enzyme_Substrates, ...
    modelParam_values.Enzyme_Products, ...
    modelParam_values.Regulators,  ...
    modelParam_values.FurtherParams];

% Sort data by reaction names:
[names_sorted, order] = sort(names_str);
values_sorted = values_array(:,order);
