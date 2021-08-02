function [F_kal] = Dynamic_factors(X,q,r,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute and determine dynamic Factors by use two-steps method
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com
% INPUTS
% X - Matrix of observable variables
% q - Number of dynamic factors
% r - Number of static factors
% p - Lags of VAR process
%
% OUTPUTS
% F_kal - Dynamic facotrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Programs')

OPTS.disp = 0;
[T,N] = size(X);
x_tot = X;
x = X;

if r < q ;
    error('Number of dynamic factors must be less than Static')
end

% Determine Static Factors by PCA
[chi, v, d, F] = PCA_estimate(x,r);
% Estimate VAR model with static facotrs
[initx, initV, A, C, R, Q] = VAR_estimate(x, chi, F, p, r, q, v);
% Parameters initialized
previous_loglik = -inf;
loglik = 0;
num_iter = 0;
LL = [];
os = size(C,1);
ss = size(A,1);
y = x';
LL = -inf;
converged = 0;
% Esimate of factors by Kalman Filter
[xitt,xittm,Ptt,Pttm,loglik_t]=Kalman_filter(initx,initV,x_tot,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=Kalman_smoother(x_tot,A,xitt,xittm,Ptt,Pttm,C,R);

F_kal =  xsmooth(1:q,:)';