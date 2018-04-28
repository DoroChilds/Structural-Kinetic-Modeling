function T = subFct_thetaMatrix(s, p, r, f, param_intervals, kinetics)
% - param_intervals: struct with fields indicating the positions of certain types of effects in the Theta matrix
% - s, p, ci, ui, r, f: Values to be assigned to the matrix
% - kinetics: type of kinetics according to which the parameters are assigned.%       
%       'enzymatic_rev':   Each pair of 2 rows is treated as forward and backward reaction so that Michaelis-Menten parameters 
%                   for row 1 are also assigned to row 2 in a modified version. (t -> t-1, -p -> 1-p)
%       other values: Michaelis-Menten parameters are jus assigned to the individual reactions for which they are declared.

Enzyme_Substrates   =   param_intervals.Enzyme_Substrates;
Enzyme_Products     =   param_intervals.Enzyme_Products;
Constant_Params =   param_intervals.Constant_Params;
Regulators      =   param_intervals.Regulators;
FurtherParams_lower     = param_intervals.FurtherParams_lower;
FurtherParams_upper     = param_intervals.FurtherParams_upper;
FurtherParams_intervals = FurtherParams_upper - FurtherParams_lower;

[n_rct, n_metab] = size(Enzyme_Substrates);

T = zeros(n_rct, n_metab);

if strcmp(kinetics, 'enzymatic_rev')
    reversible = 1;
else
    reversible = 0;
end

%% Michaelis-Menten-substrates:
[subst_MM_rows, subst_MM_cols] = find(Enzyme_Substrates);      % Michaelis-Menten substrates
n_sMM   = length(subst_MM_rows);

for iter = 1:n_sMM
    T(subst_MM_rows(iter), subst_MM_cols(iter))         = T(subst_MM_rows(iter), subst_MM_cols(iter))   + s(iter);
    if reversible
        T(subst_MM_rows(iter)+1, subst_MM_cols(iter))   = T(subst_MM_rows(iter)+1, subst_MM_cols(iter)) + s(iter) - 1;
    end
end

%% Michaelis-Menten-products:
[prods_rows, prods_cols] = find(Enzyme_Products);      % Michaelis-Menten products
n_p   = length(prods_rows);

for iter = 1:n_p
    T(prods_rows(iter), prods_cols(iter))       = T(prods_rows(iter), prods_cols(iter))     + p(iter);
    if reversible
        T(prods_rows(iter)+1, prods_cols(iter)) = T(prods_rows(iter)+1, prods_cols(iter))   + p(iter) + 1;
    end
end
    
%% Substrates of reactions following mass-action kinetics:
[subst_MWG_rows, subst_MWG_cols] = find(Constant_Params);      % Mass-action products
n_sMWG   = length(subst_MWG_rows);

for iter = 1:n_sMWG
    T(subst_MWG_rows(iter), subst_MWG_cols(iter)) = T(subst_MWG_rows(iter), subst_MWG_cols(iter)) + Constant_Params(subst_MWG_rows(iter), subst_MWG_cols(iter));
end

%% Allosteric effects:
[R_rows, R_cols] = find(Regulators);  
n_R   = length(R_rows);

for iter = 1:n_R    
    T(R_rows(iter), R_cols(iter))       = T(R_rows(iter), R_cols(iter))     + r(iter);
    if reversible
        T(R_rows(iter)+1, R_cols(iter)) = T(R_rows(iter)+1, R_cols(iter))   + r(iter);
    end
end

%% Further effects:
% no effect on the reverse reactions, even if reversible kinetics are assumed:
% (Should be used, for example, to model transport reactions where substrates and products are located in different compartments.)
[F_rows, F_cols] = find(FurtherParams_intervals);  
n_F   = length(F_rows);

for iter = 1:n_F
    T(F_rows(iter), F_cols(iter))       = T(F_rows(iter), F_cols(iter))     + f(iter);
end




