function T0 = subFct_reduce_T(T0, T1, S0, S1, L1)
% Reduce conservation reactions into the matrix Theta (saturation parameter matrix)
for dep=1:size(L1,1)
    LL = L1(dep,:) .* S0' / S1(dep); % Link matrix entry for dependent metabolite: Linear combination of all independent metabolites necessary for computing COASH-changes,
    % weighted by the ratio of concentrations ([indep. metabolites]/[dep. metabolites])
    ix = find(T1(:,dep));           % Find reactions with saturation parameters related to the dependent metabolite
    for i=1:length(ix)
        % for each reaction:
        %   Saturation parameters (reduced Theta) = Saturation parameters (reduced Theta)
        %   + Link-coeff * concentration-ratio * Saturation parameters (redundant Metabolite)
        T0(ix(i),:) =  T0(ix(i),:) + LL * T1(ix(i),dep);
    end
end