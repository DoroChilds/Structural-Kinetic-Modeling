function [label_vec, label_names] = subFct_evalEigenvalues(max_eigen)

max_real = real(max_eigen(:,1));
min_real = real(max_eigen(:,end));
label_vec = nan(size(max_real));


% detect stable states:
index_S = (max_real < -eps) & (min_real < -eps);

% detect unstable states:
index_I = (max_real > eps) & (min_real > eps);

% detect saddle points:
index_SA = (max_real > eps) & (min_real < -eps);

% models for which stability is not detectable because eigenvalues lie in the range of the machine error:
index_U = abs(max_real) < eps;

% Assign labels:
label_vec(index_I) = 0;
label_vec(index_SA) = 0;
label_vec(index_S) = 1;
label_vec(index_U) = nan;
if nargout > 1
    label_names = {'0: UNSTABLE', '1: STABLE', 'nan: Unknown'};
end
