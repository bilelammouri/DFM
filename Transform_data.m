function [x_tot,  x] = transform_data(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function transform original data to standardizied data without 
% missing value
% Author : Ammouri Bilel
% E-mail : bilel.ammouri@gmail.com
% INPUTS
% X - Matrix of observable variables (TxN)
% th - Threshold missing value
% OUTPUTS
% x_tot - Standardizied data
% x - Data without NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,N] = size(X) ;
MAT_NAN = isnan(X) ;
% Determine the maximum number of missing value in first data
nbNANmax_debut = max(sum(MAT_NAN(1:end-10,:))) ;   
% Determine the maximum number of missing value in end data
nbNANmax_fin = max(sum(MAT_NAN(end-10:end,:))) ; 
% Standardizied data
Mx= nanmean(X);  
Wx = nanstd(X); 
x_tot = (X-kron(ones(T,1),Mx))./kron(ones(T,1),Wx);
% cylinder data 
x = x_tot(nbNANmax_debut+1:end-nbNANmax_fin,:);