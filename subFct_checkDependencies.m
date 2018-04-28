function [m_dependent] = subFct_checkDependencies(options, N)
% Extract vector indicating dependent metabolites if provided by the user. 
% Also check for validity.
% If such a vector is not provided, or if the given dependencies do not fit to the number of dependent rows in N, create vector automatically.

n_met = size(N, 1);

if isfield(options, 'm_dependent')    
    m_dependent = options.m_dependent;
    % convert to column vector, if necessary:
    [d_n1, d_n2] = size(m_dependent);
    if d_n1 == 1 && d_n2 == n_met
        m_dependent = m_dependent';
    end
else
    m_dependent = zeros(n_met, 1);
end

%% Check vector indicating dependencies for validity:
% Check 1: does numer of dependent metabolites equal number of linearly dependent rows?
generic_d = 0;  % 'flag' indicating whether we need to generically generate the vector of dependencies
rank_N = rank(N);
if sum(m_dependent) ~= n_met - rank_N
    generic_d = 1;
end

% Check 2: where the correct metabolites marked as dependent?
N_reduced = N(~m_dependent,:);
if rank(N_reduced) ~= rank_N
    generic_d = 1;
end

%% Generically look for positions to remove:
if generic_d
    if options.verbose
        fprintf('Automatically detecting dependencies between metabolites.\n');
    end
    N_original = N;
    m_dependent = zeros(n_met, 1);
    count = 0;
    for row = 1:n_met
        N_temp = N;
        N_temp(row-count, :) = [];
        if rank(N_temp) == rank(N)
            m_dependent(row) = 1;
            N = N_temp;
            count = count+1;
        end
    end
    N = N_original;
end



