function [xitt,xittm,Ptt,Pttm,loglik] = Kalman_filter(initx,initV,x,A,C,R,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to use for estimate of state-space model
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com 
%INPUTS
% x(:,t) - Observation matrix at time (t) with NaN==0
% x0(:,t) - Observation matrix at time (t) with NaN  
% A - System matrix (trnasition matrix)
% C - Observation matrix (mesure matrix)
% Q - Covariance system
% R - Observation covariance
% initx - Initial state vector 
% initV - Initial state covariance
%
% OUTPUT
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
% loglik - Log-likelihood value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,N]=size(x);
r=size(A,1);
matna = isnan(x) ;
% replace NaN with random value
for i=1:1:N
	ind = find(matna(:,i)) ;
    x(ind,i) = 0 ;
end

y=x';
xittm=[initx zeros(r,T)];
xitt=zeros(r,T);
Pttm=zeros(r,r,T);
Pttm(:,:,1)=initV;
Ptt=zeros(r,r,T);
variance = diag(R); 

for j=1:T
    % Preprocessing missing value
    i_d = find(matna(j,:)); 
    variance_d = variance; 
    variance_d(i_d) = 10^10; 
    R=diag(variance_d);
    L=inv(C*Pttm(:,:,j)*C'+R);

    
    %% Klman Filtre steps
    % Steps 1 and 5
    xitt(:,j)=xittm(:,j)+Pttm(:,:,j)*C'*L*(y(:,j)-C*xittm(:,j));
    % Steps 2
    Ptt(:,:,j)=Pttm(:,:,j)-Pttm(:,:,j)*C'*L*C*Pttm(:,:,j); 
    % Steps 3
    xittm(:,j+1)=A*xitt(:,j);  
    % Steps 4
    Pttm(:,:,j+1)=A*Ptt(:,:,j)*A'+Q; 
    
    e = y(:,j) - C*xittm(:,j);
    n = length(e);
    ss = length(A);
    d = size(e,1);
    S = C*Pttm(:,:,j)*C' + R;
    GG = C'*diag(1./diag(R))*C;
    Sinv = inv(S);

    detS = prod(diag(R))*det(eye(ss)+Pttm(:,:,j)*GG);
    denom = (2*pi)^(d/2)*sqrt(abs(detS));
    mahal = sum(e'*Sinv*e,2);
    logl(j) = -0.5*mahal - log(denom);
    
end

loglik=sum(logl);