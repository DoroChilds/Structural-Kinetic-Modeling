function [RegulatorMatrix, positions, values] = skm_assignRandomRegulators(networkDimensions, intervalLimits, kinetics)
% Creates a matrix with regulatory parameters at random network positions that can be 
% assigned to the 'Regulators'-field of a 'paramInterval'-struct

if nargin < 3
    kinetics = 'enzymatic_irrev';
end

number_of_regulators = length(intervalLimits);
RegulatorMatrix      = zeros(fliplr(networkDimensions));

if ~ strcmp(kinetics, 'enzymatic_rev')
    number_of_matrix_elements   = numel(RegulatorMatrix);
    shuffled_matrix_positions   = randperm(number_of_matrix_elements);
    regulator_positions         = shuffled_matrix_positions(1:number_of_regulators);
    RegulatorMatrix(regulator_positions) = intervalLimits;
else
    RegulatorMatrix_forwardRct_only = RegulatorMatrix(1:2:end,:);
    number_of_matrix_elements       = numel(RegulatorMatrix_forwardRct_only);
    shuffled_matrix_positions       = randperm(number_of_matrix_elements);
    regulator_positions             = shuffled_matrix_positions(1:number_of_regulators);
    RegulatorMatrix_forwardRct_only(regulator_positions) = intervalLimits;
    RegulatorMatrix(1:2:end,:)      = RegulatorMatrix_forwardRct_only;
end

[regulator_rows, regulator_cols, values] = find(RegulatorMatrix);
positions = [regulator_rows, regulator_cols];