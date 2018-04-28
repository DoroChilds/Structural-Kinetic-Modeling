function [eigen_matrix, params_used, param_names] = subFct_performIteration(N, dependent_metabolites, S0, v0, param_intervals, num_models, options, m_names, r_names)

%% Some preliminary computations
[n_metab, n_rct] = size(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine number of parameters to sample: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[param_info, param_index] = subFct_thetaNumbers(param_intervals);
n_S     = param_info.n_S;
n_P     = param_info.n_P;
n_R     = param_info.n_R;
n_F     = param_info.n_F;
r_limits        = param_info.r_limits;
f_intervals     = param_info.f_intervals;
f_lowerbounds   = param_info.f_lowerbounds;
numParams = n_S+n_P+n_R+n_F;
if options.verbose
    fprintf('Number of sampled model parameters: %g\n', numParams)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly sample parameters: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options, 'simType')
    simType = options.simType;
else
    simType = 'MonteCarlo';
end
if strcmp(simType, 'MonteCarlo')
    params_s    = rand(num_models,n_S);
    params_p    = -rand(num_models,n_P);
    params_r    = rand(num_models,n_R)  .*  repmat(r_limits,num_models,1);
    params_f    = rand(num_models,n_F)  .*  repmat(f_intervals,num_models,1) + repmat(f_lowerbounds, num_models, 1);
elseif strcmp(simType, 'gridSearch')
    % Preform grid search:
    points_per_param = num_models;
    grid_vec = linspace(0+eps,1-eps,points_per_param)';
    gridMatrix = grid_vec;
    for p = 1:numParams-1
        current_vec1 = repmat(gridMatrix(:,1),1,points_per_param)';
        current_vec2 = current_vec1(:);
        gridMatrix_rep = repmat(gridMatrix, points_per_param, 1);
        gridMatrix = [current_vec2, gridMatrix_rep];
    end
    num_models = size(gridMatrix,1);
    params_s = gridMatrix(:,1:n_S); gridMatrix(:,1:n_S)=[];
    params_p = gridMatrix(:,1:n_P); gridMatrix(:,1:n_P)=[];
    params_r = gridMatrix(:,1:n_R); gridMatrix(:,1:n_R)=[];
    params_f = gridMatrix(:,1:n_F); gridMatrix(:,1:n_F)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve indices for the different compound and reaction names %
% ( will be passed on to the optional function fct_modify_params() ):  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_indices = subFct_names_to_indices(m_names, r_names);

%% Cope with dependent metabolites
% Compute Link matrix:
[N_r N0 N1 dependent_metabolites] = subFct_rearrangeN(N, dependent_metabolites);
[L L0 L1]   = subFct_linkMatrix(N0, N1);

% Sort metabolites so that they fit to link matrix:
[S_r S0_0 S0_1] = subFct_rearrangeN(S0, dependent_metabolites);
% Test: sum(sum(abs(L * N0 - N_r)))

%% Create matrix of eigenvalues
if options.return_all_eigvals
    eigen_matrix = nan(num_models, size(N0,1));
else
    eigen_matrix = nan(num_models, 1);
end
    
%% Start iteration
if options.verbose
%     fprintf('\n');
end
for iter=1:num_models
    if options.verbose && mod(iter,10^3)==0
        fprintf('Iterations: %g\n', iter);
    end

    s_current_params = params_s(iter,:);
    p_current_params = params_p(iter,:);
    r_current_params = params_r(iter,:);
    f_current_params = params_f(iter,:);

    T = subFct_thetaMatrix(s_current_params, p_current_params, r_current_params, f_current_params, param_intervals, options.kinetics);
    T = feval(options.fct_modify_params, T, model_indices, options);
    
    [T_r T0 T1 dependent_metabolites] = subFct_rearrangeN(T', dependent_metabolites);
    T_r = T_r'; 
    T0 = T0'; 
    T1 = T1';
    
    % Normalize:
    Lambda = N .* repmat(v0', n_metab, 1) ./ repmat(S0, 1, n_rct);
    [Lambda_r Lambda0 Lambda1 dependent_metabolites] = subFct_rearrangeN(Lambda, dependent_metabolites);
    
    % Reduce conservation reactions:
    T0 = subFct_reduce_T(T0, T1, S0_0, S0_1, L1);

    % Compute Jacobian:
    J = Lambda0 * T0;
    
    % Compute eigenvalues:
    if options.return_all_eigvals
        %         tic
        eigen = eig(J);        
        % Sort according to real parts (maximum real part is first) and store in output matrix:
        [real_values_sorted, sort_index] = sort(real(eigen), 'descend');
        eigen_matrix(iter,:) = eigen(sort_index);
        %         toc
    else
        %         tic
        eigen_matrix(iter) = eigs(J, 1, 'lr');
        %         toc
    end
    
end
if options.verbose
    %     fprintf('\n');
end

%% Prepare output

% Assign names to the model parameters:
param_names.Enzyme_Substrates   = [num2cell(1: size(param_index.s,1))', r_names(param_index.s(:,1))', m_names(param_index.s(:,2))', repmat({'Substrate'}, size(param_index.s,1), 1)];
param_names.Enzyme_Products     = [num2cell(1: size(param_index.p,1))', r_names(param_index.p(:,1))', m_names(param_index.p(:,2))', repmat({'Product'},   size(param_index.p,1), 1)];
param_names.Regulators      = [num2cell(1: size(param_index.r,1))', r_names(param_index.r(:,1))', m_names(param_index.r(:,2))', repmat({'Regulator'}, size(param_index.r,1), 1)];
param_names.FurtherParams   = [num2cell(1: size(param_index.f,1))', r_names(param_index.f(:,1))', m_names(param_index.f(:,2))', repmat({'Further_Parameter'}, size(param_index.f,1), 1)];

% Collect sampled parameters:
params_used.Enzyme_Substrates   = params_s;
params_used.Enzyme_Products     = params_p;
params_used.Regulators      = params_r;
params_used.FurtherParams   = params_f;

% Sort parameters according to their reaction names:
param_types = fieldnames(param_names);
for iter = 1:length(param_types)
    temp_type = param_types{iter};
    % Sort names:
    temp_names = param_names.(temp_type);
    [names_sorted_temp, order] = sort(temp_names(:,2));  % hier weiter mit debuggen!
    param_names.(temp_type) = temp_names(order,2:end);
    % Sort parameters:
    temp_values = params_used.(temp_type);
    params_used.(temp_type) = temp_values(:, order);
end







