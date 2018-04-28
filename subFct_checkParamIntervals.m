function param_intervals = subFct_checkParamIntervals(user_input, N, kinetics)

[n_met, n_rct] = size(N);
param_intervals = user_input;

if ~isfield(param_intervals, 'Enzyme_Substrates')
    
    if strcmp(kinetics, 'enzymatic_irrev')
        param_intervals.Enzyme_Substrates = double(N<0)';

    elseif strcmp(kinetics, 'enzymatic_rev')
        % Split matrix into forward- and backward columns:
        N_forward = N(:,1:2:end-1);
        N_backward = N(:,2:2:end);
        % Find substrates of forward reactions for which the backward reaction has a complementary stoichiometric coefficient:
        N_forward_MM_substrates = double(N_forward == -N_backward);
        N_forward_MM_substrates(N_forward>=0) = 0;
        % Mark positions of Michaelis-Menten-substrates in a matrix with the original dimensions:
        N_complete_MM_substrate = zeros(n_met, n_rct);
        N_complete_MM_substrate(:,1:2:end-1) = N_forward_MM_substrates;
        % Save positions of Michaelis-Menten-substrate parameters:
        param_intervals.Enzyme_Substrates = N_complete_MM_substrate';
        
    else
        param_intervals.Enzyme_Substrates = zeros(n_rct, n_met);
    end
end
if ~isfield(param_intervals, 'Enzyme_Products')
    if strcmp(kinetics, 'enzymatic_rev')
        % Split matrix into forward- and backward columns:
        N_forward = N(:,1:2:end-1);
        N_backward = N(:,2:2:end);
        % Find products of forward reactions for which the backward reaction has a complementary stoichiometric coefficient:
        N_forward_MM_products = double(N_forward == -N_backward);
        N_forward_MM_products(N_forward<=0) = 0;
        % Mark positions of Michaelis-Menten-substrates in a matrix with the original dimensions:
        N_complete_MM_product = zeros(n_met, n_rct);
        N_complete_MM_product(:,1:2:end-1) = N_forward_MM_products;
        % Save positions of Michaelis-Menten-substrate parameters:
        param_intervals.Enzyme_Products = N_complete_MM_product';
        
    else
        param_intervals.Enzyme_Products = zeros(n_rct, n_met);
    end
end

if ~isfield(param_intervals, 'Constant_Params')
    if  strcmp(kinetics, 'massAction')
        param_intervals.Constant_Params = N'<0;       % check later: 1/N?
    else
        param_intervals.Constant_Params = zeros(n_rct, n_met);
    end
end
if ~isfield(param_intervals, 'Regulators')
    param_intervals.Regulators = zeros(n_rct, n_met);
end
if ~isfield(param_intervals, 'FurtherParams_lower')
    param_intervals.FurtherParams_lower = zeros(n_rct, n_met);
end
if ~isfield(param_intervals, 'FurtherParams_upper')
    param_intervals.FurtherParams_upper = zeros(n_rct, n_met);
end

