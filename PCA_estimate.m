function [chi, v, d, F] = PCA_estimate(X,r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute and determine the Static Factors
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com
% INPUTS
% X - Matrix of observable variables (TxN)
% r - Number of static factors
%
% OUTPUTS
% F - Static facotrs
% chi - Commun component gamma_chap * Ft   (T x N)
% v - Eigenvector matrix
% d - Enginevalue matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPTS.disp = 0;
[x] = Transform_data(X);
[v, d] = eigs(cov(x),r,'lm',OPTS); 

chi = x*v*v';
d = eye(r);
F=x*v; 