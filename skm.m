function [eigenvalues, stability, modelParam_values, modelParam_names, used_options, used_paramIntervals] = skm(N, S0, v0, numModels, options, paramIntervals)
% 
% [eigenvalues, stability, modelParam_values, modelParam_names, used_options, used_paramIntervals] = ...
% skm(N, S0, v0, numModels, options, paramIntervals)
%
% The main function of the SKM-toolbox. It creates the SK-models, performs Monte-Carlo 
% simulations and plots the resulting distributions of eigenvalues and model parameters.

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set default argumenst, if necessary %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6
    paramIntervals = struct();
    if nargin < 5
        options = struct();
        if nargin < 4
            numModels = 10^4;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check options for consistency: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = subFct_checkOptions(options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve / create model-specific data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------- %
% Model components: %
% ----------------- %
% Model size:
[n_met, n_rct] = size(N);

% Names of metabolites and reactions: use user-defined values or otherwise create them generically.
[m_names, r_names] = subFct_checkNames(options, n_met, n_rct);

% Indices of the dependent metabolites: use user-defined vector, or create it generically.
options.m_dependent = subFct_checkDependencies(options, N);
used_options        = options;


% ----------------- %
% Model parameters: %
% ----------------- %
% use user-defined matrices, or create them generically according to the desired kinetic rate laws:
paramIntervals = subFct_checkParamIntervals(paramIntervals, N, options.kinetics);
used_paramIntervals = paramIntervals;

% ------------------ %
% Steady state data: %
% ------------------ %
% convert to column vectors, if necessary:
[r_n1, r_n2] = size(v0);
[m_n1, m_n2] = size(S0);
if r_n1 == 1 && r_n2 == n_rct
    v0 = v0';
end
if m_n1 == 1 && m_n2 == n_met
    S0 = S0';
end

% ---------------------------- %
% Check model for consistency: %
% ---------------------------- %
% Check constistency:
consistent = subFct_checkModelConsistency(N, v0, S0, paramIntervals, options.verbose);
% If errors exist, create 'dummy' return values and terminate this function:
if ~consistent
    eigenvalues = []; modelParam_names = []; modelParam_values = []; stability = [];
    fprintf('Terminating function due to errors!\n');
    return
end


% -------------------------------------------- %
% Remove unnecessary metabolites and reactions %
% -------------------------------------------- %
% Remove components desired by the user:
[N, S0, m_dependent, m_names, paramIntervals] = subFct_removeMetabolites(options.m_exclude, options.m_include, N, S0, options.m_dependent, paramIntervals, m_names, options.verbose);
[N, v0, r_names, paramIntervals]              = subFct_removeReactions  (options.r_exclude, options.r_include, N, v0, paramIntervals, r_names, options.verbose);

% Remove components through automatic checks (optional):
[N, S0, v0, m_dependent, paramIntervals, m_names, r_names] = subFct_removeZeroElements(N, S0, v0, m_dependent, paramIntervals, ...
    m_names, r_names, ...
    options.rm_zero_met, options.rm_zero_rct, options.rm_zero_N, options.verbose, options.kinetics);

% Test whether steady state assumption is still fulfilled with reduced N and v:
dSdt    = N * v0;
max_dSdt = max(abs(dSdt));
if options.verbose
    fprintf('Maximum concentration change in the steady state: max(abs(dS/dt)) = %2.2e\n', max_dSdt);
end


%%%%%%%%%%%%%%%%%%%
%% Start SK-model %
%%%%%%%%%%%%%%%%%%%
% ------------------------------- %
% Build desired number of models: %
% ------------------------------- %
[eigenvalues, modelParam_values, modelParam_names] = subFct_performIteration(N, m_dependent, S0, v0, paramIntervals, numModels, options, m_names, r_names);
% numModels = size(eigenvalues, 1);

% --------------------------------------------------------- %
% Balance number of stable and unstable classes (optional): %
% --------------------------------------------------------- %
if options.balance_stability
    [L_STABLE, L_NAMES_STABLE] = subFct_evalEigenvalues(eigenvalues);
    models_per_class = min([sum(L_STABLE==1), sum(L_STABLE==0)]);
    
    % Repeat the model production until enough stable and unstable models exist:
    while models_per_class < numModels/2
        fprintf('\nRepeating sampling procedure until equal numbers\nof stable and unstable models are reached!\n\n')
        [EIGEN_new, S_params_new, param_index_new] = subFct_performIteration(N, m_dependent, S0, v0, paramIntervals, numModels, options, m_names, r_names);
        modelParam_values.Enzyme_Substrates = [modelParam_values.Enzyme_Substrates; S_params_new.Enzyme_Substrates];
        modelParam_values.Enzyme_Products = [modelParam_values.Enzyme_Products; S_params_new.Enzyme_Products];
        modelParam_values.Regulators = [modelParam_values.Regulators; S_params_new.Regulators];
        modelParam_values.FurtherParams = [modelParam_values.FurtherParams; S_params_new.FurtherParams];
        eigenvalues = [eigenvalues; EIGEN_new];
        [L_STABLE, L_NAMES_STABLE] = subFct_evalEigenvalues(eigenvalues);
        models_per_class = min([sum(L_STABLE==1), sum(L_STABLE==0)]);
    end
    
    % Produce balanced dataset:
    i_C1 = find(L_STABLE==1);
    i_C1 = i_C1(1:numModels/2);
    i_C0 = find(L_STABLE==0);
    i_C0 = i_C0(1:numModels/2);
    modelParam_values.Enzyme_Substrates = [modelParam_values.Enzyme_Substrates(i_C1,:); modelParam_values.Enzyme_Substrates(i_C0,:)];
    modelParam_values.Enzyme_Products = [modelParam_values.Enzyme_Products(i_C1,:); modelParam_values.Enzyme_Products(i_C0,:)];
    modelParam_values.Regulators = [modelParam_values.Regulators(i_C1,:); modelParam_values.Regulators(i_C0,:)];
    modelParam_values.FurtherParams = [modelParam_values.FurtherParams(i_C1,:); modelParam_values.FurtherParams(i_C0,:)];
    eigenvalues = [eigenvalues(i_C1,:); eigenvalues(i_C0,:)];
    
    % Mix stable and unstable models with random order:
    i_rand = randperm(numModels);
    modelParam_values.Enzyme_Substrates = modelParam_values.Enzyme_Substrates(i_rand,:);
    modelParam_values.Enzyme_Products = modelParam_values.Enzyme_Products(i_rand,:);
    modelParam_values.Regulators = modelParam_values.Regulators(i_rand,:);
    modelParam_values.FurtherParams = modelParam_values.FurtherParams(i_rand,:);
    eigenvalues = eigenvalues(i_rand,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare return values %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign labels according to stability:
[stability, L_NAMES_STABLE] = subFct_evalEigenvalues(eigenvalues);
[pmean, psd] = skm_countStableModels(stability);
[Umean, Usd] = skm_countStableModels(stability==0);
[Nmean, Nsd] = skm_countStableModels(isnan(stability));

%%%%%%%%%%%%%%%%%%
%% Print results %
%%%%%%%%%%%%%%%%%%
if options.verbose
    S_index = stability == 1;
    I_index = stability == 0;
    U_index = isnan(stability);
    
    fprintf('Results:\n')
    fprintf('Percentage of stable models (max_eig < %1.2e): %2.2f +- %2.2f %%\n', eps, pmean, psd)
    fprintf('Percentage of unstable models (max_eig > %1.2e): %2.2f +- %2.2f %%\n' ,   eps, Umean, Usd)
    fprintf('Percentage of unclear models (|max_eig| <= %1.2e): %2.2f +- %2.2f %%\n' ,   eps, Nmean, Nsd)
    
    time_stop = toc;
    fprintf('Elapsed time: %5.2f seconds', time_stop);    
    fprintf('\n\n')
end

%%%%%%%%%%%%%%%%%
%% Plot results %
%%%%%%%%%%%%%%%%%
if options.plot
    h_vec = skm_plotParams(modelParam_values, eigenvalues, stability, modelParam_names);
end
