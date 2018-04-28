function consistent = subFct_checkModelConsistency(N, v0, S0, param_intervals, verbose)
% This subroutine checks the given model components (stoichiometry, steady state data, matrices with parameter intervals) for consistency.

consistent = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if steady state data are given as vectors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_n1, r_n2] = size(v0);
[m_n1, m_n2] = size(S0);
if min([r_n1, r_n2]) > 1
    consistent = 0;
    if verbose
        fprintf('Error: Steady state reaction rates are not given in a vector, but in a matrix of size %gx%g. Please remove extra columns or rows!\n', r_n1, r_n2);
    end
end
if min([m_n1, m_n2]) > 1
    consistent = 0;
    if verbose
        fprintf('Error: Steady state concentrations are not given in a vector, but in a matrix of size %gx%g. Please remove extra columns or rows!\n', m_n1, m_n2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% Compare dimensions: %
%%%%%%%%%%%%%%%%%%%%%%%
[n_met, n_rct] = size(N);

if length(v0) ~= n_rct
    consistent = 0;
    if verbose
        fprintf('Error: Number of steady state reaction rates (%g) does not correspond to column number of the stoichiometric matrix (%g). Please check again!\n', length(v0), n_rct);
    end
end
if length(S0) ~= n_met
    consistent = 0;
    if verbose
        fprintf('Error: Number of steady state concentrations (%g) does not correspond to row number of the stoichiometric matrix (%g). Please check again!\n', length(S0), n_met);
    end
end

param_types = fieldnames(param_intervals);
for i = 1:length(param_types)
    current_param_matrix = param_intervals.(param_types{i}); 
    [n1, n2] = size(current_param_matrix);
    if n1~=n_rct || n2~=n_met
        consistent = 0;
        if verbose
            fprintf('Error: Parameter-interval matrix "%s" has wrong size (%gx%g). Please make sure that it corresponds to the size of the transposed stoichiometric matrix (%gx%g)!\n', param_types{i}, n1, n2, n_rct, n_met);
        end
    end
end

