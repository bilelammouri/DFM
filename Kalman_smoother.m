function [xitT,PtT,PtTm]=Kalman_smoother(x,A,xitt,xittm,Ptt,Pttm,C,R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to use for lissage of state-space model
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com
% INPUTS
% y(:,t) - Observation at time t
% A - System matrix
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]  (taille r x r x N)
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]  (taille r x r x N)
% C - Observation matrix 
% R - Observation covariance
%
% OUTPUT:
% xitT = E[X(:,t) | y(:,1:T)]
% PtT = Cov[X(:,t) | y(:,1:T)]
% PtTm = Cov[X(:,t+1),X(:,t) | y(:,1:T)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,N]=size(x);
r=size(A,1);
Pttm=Pttm(:,:,1:end-1);
xittm=xittm(:,1:end-1);
J=zeros(r,r,T);
variance = diag(R);
matna = isnan(x) ;

% Replace missing value with random value
for i=1:1:N
	ind = find(matna(:,i)) ;
    x(ind,i) = 0;
end

for i=1:T-1
    J(:,:,i)=Ptt(:,:,i)*A'*inv(Pttm(:,:,i+1));
end

for i=1:T
    % Preprocessing missing value
    L(:,:,i)=inv(C*Pttm(:,:,i)*C'+R);
    K(:,:,i)=Pttm(:,:,i)*C'*L(:,:,i);
end

xitT=[zeros(r,T-1)  xitt(:,T)];
PtT=zeros(r,r,T);
PtTm=zeros(r,r,T);
PtT(:,:,T)=Ptt(:,:,T);
PtTm(:,:,T)=(eye(r)-K(:,:,T)*C)*A*Ptt(:,:,T-1);
%% Lissage function
for j =1:T-1
    % Equation 1 : expectation of factors
    xitT(:,T-j)= xitt(:,T-j)+J(:,:,T-j)*(xitT(:,T+1-j)-xittm(:,T+1-j));
    % Equation 2 : variance of factors Z_(t+1)/T
    PtT(:,:,T-j)=Ptt(:,:,T-j)+J(:,:,T-j)*(PtT(:,:,T+1-j)-Pttm(:,:,T+1-j))*J(:,:,T-j)';   
end

for j =1:T-2
    % Equation 2 : variance of factors Z_(t)/T
    PtTm(:,:,T-j)=Ptt(:,:,T-j)*J(:,:,T-j-1)'+J(:,:,T-j)*(PtTm(:,:,T-j+1)-A*Ptt(:,:,T-j))*J(:,:,T-j-1)'; 
end