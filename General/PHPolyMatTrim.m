function [R2,rho] = PHPolyMatTrim(R1,gamma)
%[R2,rho] = PHPolyMatTrim(R1,gamma)
% 
%  R2=PHPolyMatTrim(R1) trims the lag dimension of a parahermitian matrix R1 by
%  symmetrically removing the outer matrix coefficients. The function will 
%  remove the matrix coefficient with the smallest Frobenius
%  norm, such that by the total energy of the trimmed parts does not exceed  
%  one permille of the total energy in R1.
% 
%  R2=PHPolyMatTrim(R1,gamma) trims such that the ratio between the removed 
%  and total energy of R1 is less than gamma, with 0<=gamma<1.
%
%  [R2,rho]=PUPolyMatTrim(R1) or [R2,rho]=PUPolyMatTrim(R1,r) additionally 
%  returns the ratio of the actual suppressed energy in R2.
%  
%  This function is related to PUPolyMatTrim() and the 'trim' operation 
%  described in [1,2].
%
%  Input parameters:
%     R1      parahermitian matrix
%     gamma   maximum ratio between the maximum energy removed at outer lags 
%             and the total energy 
%             default: 1/1000
%  
%  Output parameters:
%     R2      trimmed parahermitian matrix
%     rho     ratio of removed energy
%
%  References:
%  [1] J.G. McWhirter, P.D. Baxter, T. Cooper, S. Redif, and J. Foster, "An EVD 
%      Algorithm for Para-Hermitian Polynomial Matrices," IEEE Transactions on 
%      Signal Processing, vol. 55, no. 5, pp. 2158-2169, May 2007.
%  [2] C.H. Ta and S. Weiss, "Shortening the Order of Paraunitary Matrices in  
%      SBR2 Algorithm", 6th International Conference on Information, Communi-
%      cations & Signal Processing, Singapore, pp. 1-5, Dec. 2007.

% S. Weiss and J. Corr, University of Strathclyde, 14/11/2014

% check input parameters
if nargin<2,
   gamma = 0.001;
end;

% check limit
Norm1 = PolyMatNorm(R1);
E = gamma*Norm1;
[M,~,L] = size(R1);

% trim from sequentially from front or back, where ever the matrix with the 
% smallest norm can be found, until the desired amount of energy has been 
% removed
Etrim = 0;
IndexRemove = 1;
Efront = norm(R1(:,:,IndexRemove),'fro').^2;
while ((Etrim + 2*Efront) <= E),
   Etrim = Etrim + 2*Efront;
   IndexRemove = IndexRemove+1;
   Efront = norm(R1(:,:,IndexRemove),'fro').^2;
end;    

% output
R2 = R1(:,:,IndexRemove:L+1-IndexRemove);
Etrim = Etrim - 2*Efront;
rho = Etrim/Norm1;
