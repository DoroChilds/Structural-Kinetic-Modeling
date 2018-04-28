function [p_values_sorted, names_sorted] = skm_pairwiseTests(modelParam_values, stability, modelParam_names)
% [p_values_sorted, names_sorted] = skm_pairwiseTests(modelParam_values, stability, modelParam_names)
%
% Pairwise comparison between the model parameters responsible for stable and unstable 
% models. Uses the Kolmogorov-Smirnov test.


% Concatenate all parameters to one matrix and their names to one cell:
[params, names] = subFct_concatenateParams(modelParam_values, modelParam_names);
[n_models, n_params] = size(params);

% Distinguish between stable and unstable models:
S_index = stability==1;     % Stable model positions (S = Stability)
I_index = stability==0;     % Unstable models positions (I = Instability)
params_S = params(S_index,:);   % Parameters producing stable models
params_I = params(I_index,:);   % Parameters producing unstable models

% Original version used by Grimbs et al. (2007): Test unstable models against background distribution.
% This is not done here, because the production of balanced classes (if activated) can lead to 'biased' non-uniform distributions in 'params'
% params_S = params;

% If both classes exist...
p_values = zeros(n_params, 1);
for i = 1:n_params
    currentPar_S = params_S(:,i);
    currentPar_I = params_I(:,i);
    % Kolmogorov-Smirnov Test:
    [h_temp, p_temp] = kstest2(currentPar_S, currentPar_I);
    p_values(i,1) = p_temp;
end


[p_values_sorted, order] = sort(p_values, 'ascend');
names_sorted = names(order);
