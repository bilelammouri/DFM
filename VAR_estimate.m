function [initx, initV, A, C, R, Q] = VAR_estimate(x, chi, F, p, r, q, v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function use static factors to modelised via VAR model and estiamte 
% the parameters model
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com
% INPUTS
% x - Matrix of observable variables (TxN)
% q - Nomber of dynamic factors
% r - Nomber of static factors
% p - lags of VAR process (commun factors)
% chi - Commun component
% F - Static facotrs
% v - Eigenvector matrix
% 
% OUTPUTS
% initx  - Initial state (column) vector   
% initV -  Initial state covariance 
% A - System matrix : matrice de transition
% C - Observation matrix : matrice de mesure
% R - Observation covariance : variance du vecteurs des innovations
% Q - System covariance : variance du vecteurs des erreurs de mesure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPTS.disp = 0;
[T,N] = size(x);
nlag = p-1;  

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag+1));
A = [A_temp';I(1:end-r,1:end)];
Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

if p > 0    
    z = F;
    Z = [];
    for kk = 1:p
        Z = [Z z(p-kk+1:end-kk,:)]; 
    end;
    z = z(p+1:end,:);
    %% Use VAR model to estimate factors of transition matrix
    % chi_(t) = A*chi_(t-1) + e_(t);
    A_temp = inv(Z'*Z)*Z'*z;
    A(1:r,1:r*p) = A_temp';
    e = z  - Z*A_temp;
    % Var-cov matrix of residual VAR model
    H = cov(e);   
    if r > q 
        % Use PCA to extract q Eigenvector and Eigenvalue
        [P, M] = eigs(H,q,'lm',OPTS);    
        P = P*diag(sign(P(1,:)));
        % The common shocks
        u_orth = e*P*(M^-.5);       
        e_pc = e*P*P';
        % Chocs Variance of VAR model
        Q(1:r,1:r) = P*M*P';        
    else        
        Q(1:r,1:r) = H;
    end;
end;

R = diag(diag(cov(x-chi)));

z = F;
Z = [];
for kk = 0:nlag
    Z = [Z z(nlag-kk+1:end-kk,:)];
end;

initx = Z(1,:)';                                                            
initV = reshape(pinv(eye((r*(nlag+1))^2)-kron(A,A))*Q(:),r*(nlag+1),r*(nlag+1));eye(r*(nlag+1));

C = [v zeros(N,r*(nlag))];
