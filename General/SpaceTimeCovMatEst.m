function Rhat = SpaceTimeCovMatEst(X,MaxLag);
%SpaceTimeCovMatEst(X); 
%
%  Rhat = SpaceTimeCovMatEst(X,MaxLag) estimates the space-time covariance 
%  matrix over the lag range (-MagLag ... +MaxLag) based on L snapshots of 
%  M-array data. The estimate Rhat is unbias and has a variance as described 
%  in [1].
%
%  Input parameters
%     X           MxL  data matrix
%     MaxLag      described desired lag range of Rhat
% 
%  Output parameter
%     Rhat        estimated space-time covariance matrix
%
%  Reference:
%
%  [1] C. Delaosa, J. Pestana, and S. Weiss: "..." ICASSP'19.

% S. Weiss, UoS, 14/10/18

%------------------------------------------------------------------------------
%  parameters
%------------------------------------------------------------------------------
[M,L] = size(X);
Rhat = zeros(M,M,2*MaxLag+1);

%------------------------------------------------------------------------------
%  evaluation
%------------------------------------------------------------------------------
Rhat(:,:,MaxLag+1) = X*X'/L; 
for tau = 1:MaxLag,
  Rhat(:,:,MaxLag+1+tau) = X(:,1+tau:L)*X(:,1:L-tau)'./(L-tau);
  Rhat(:,:,MaxLag+1-tau) = Rhat(:,:,MaxLag+1+tau)';
end;

