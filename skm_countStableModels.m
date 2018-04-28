function [stablePercent_mean, stablePercent_sd, stableModels_mean, stableModels_sd, stablePercent_perRun, stableModels_perRun] = skm_countStableModels(stability, k)
% Quantitative evaluation of the number of stable and unstable models. Computes
% absolute numbers and percentages. Returns estimates and standard deviations.

if nargin < 2
    k = 10; % number of runs
end

if length(stability)>1
    %% Count mean and standard deviation for the number of stable and unstable steady states
    % In order to obtain standard deviations, the MATLAB function 'crossval' is used to
    % partition the data into 10 subsets and to compute seperate counts per subset.
    
    % Function to compute the number of stable models:
    countstable = @(x, y) sum(y==1);
    
    % Evalute results per subset:    
    % 1. absolute numbers of stable models
    stableModels_perRun = crossval(countstable, stability, 'kfold', k);
    stableModels_mean  = mean(stableModels_perRun);
    stableModels_sd    = std(stableModels_perRun);
    % 2. relative percentage proportions of stable models
    stablePercent_perRun = stableModels_perRun/(length(stability)/k) *100;
    stablePercent_mean  = mean(stablePercent_perRun);
    stablePercent_sd    = std(stablePercent_perRun);
    
else
    stableModels_mean  = stability==1;
    stableModels_sd    = 0;
end



if length(stability)>1
else
    stablePercent_mean = (stability==1)*100;
    stablePercent_sd   = 0;
end

