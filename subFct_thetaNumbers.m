function [param_info, param_index] = subFct_thetaNumbers(param_intervals)

Enzyme_Substrates   =   param_intervals.Enzyme_Substrates;
Enzyme_Products     =   param_intervals.Enzyme_Products;
Regulators      =   param_intervals.Regulators;
FurtherParams_lower     = param_intervals.FurtherParams_lower;
FurtherParams_upper     = param_intervals.FurtherParams_upper;
FurtherParams_intervals = FurtherParams_upper - FurtherParams_lower;

%% Find indices of saturation parameters:
index_S     = find(Enzyme_Substrates)';     % Michaelis-Menten substrates
index_P     = find(Enzyme_Products)';       % Michaelis-Menten products
index_R     = find(Regulators)';        % Activators
index_F     = find(FurtherParams_intervals)'; % Further parameters


%% Define numbers of parameters of each type:
n_S     = length(index_S);
n_P     = length(index_P);
n_R     = length(index_R);
n_F     = length(index_F);

%% User-defined interval boundaries:
r_limits        = Regulators(index_R);
f_intervals     = FurtherParams_intervals(index_F);
f_lowerbounds   = FurtherParams_lower(index_F);

%% Assign row and column indices:
[S_rows, S_cols]    = find(Enzyme_Substrates);
[P_rows, P_cols]    = find(Enzyme_Products);
[R_rows, R_cols]    = find(Regulators);
[F_rows, F_cols]    = find(FurtherParams_intervals);

param_index.s = [S_rows, S_cols];
param_index.p = [P_rows, P_cols];
param_index.r = [R_rows, R_cols];
param_index.f = [F_rows, F_cols];

param_info.n_S = n_S;
param_info.n_P = n_P;
param_info.n_R = n_R;
param_info.n_F = n_F;
param_info.r_limits         = r_limits;
param_info.f_intervals      = f_intervals;
param_info.f_lowerbounds    = f_lowerbounds;
