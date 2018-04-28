function [N_r N0 N1 I] = subFct_rearrangeN(N, I)
index_independent   = find(I==0);
index_dependent     = find(I~=0);

N0 = N(index_independent,:);
N1 = N(index_dependent, :);

N_r = [N0; N1];