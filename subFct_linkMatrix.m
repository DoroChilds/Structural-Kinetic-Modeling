function [L L0 L1] = subFct_linkMatrix(N0, N1)

L1 = round(N1 * pinv(N0));
L0 = eye(size(N0,1));

L = [L0; L1];