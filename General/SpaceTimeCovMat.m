function [R,tau] = SpaceTimeCovMat(X,MaxLag);
%[R,tau] = SpaceTimeCovMat(X,MaxLag);
%
%   R=SpaceTimeCovMat(X,MaxLag) estimates the polynomial covariance matrix
%   R[tau] from M-channel data x[n] such that
%       R[tau] = E{ x[n] x^H[n-tau] } 
%   Each of the M rows of X represents a time series, and each column 
%   represents the M-element vector x[n] at a specific sample index n. The 
%   time series are assumed to have zero mean such that R is a covariance
%   rather than a correlation matrix.
%
%   The evaluation of the covariance is restricted to lag values within the 
%   interval [-Maxlag;+MaxLag]. The output can also be interpreted as a 
%   polynomial or cross spectral density matrix
%      R(z) = R_{-MaxLag} z^MaxLag + ... + R_{-1} z + R0 +
%                  + R[1]z^{-1} + ... + R[MaxLag]z^{-MaxLag}
%   whereby the returned format is
%      R(:,:,1) = R_{-MagLag};
%      ...
%      R(:,:,MaxLag) = R_{-1};
%      R(:,:,MaxLag+1) = R_{0};
%      R(:,:,MaxLag+2) = R_{1};
%      ...
%      R(:,:,2*MaxLag+1) = R_{Maxlag}
%   The space-time covariance matrix is parahermitian, such that
%      R(:,:,1) = R(:,:,2*MaxLag+1)';
%      ...
%      R(:,:,MaxLag) = R(:,:,MaxLag+2)';
%   holds.
%
%   [R,tau]=SpaceTimeCovMat(X,MaxLag) additionally returns the lag para-
%   meters tau=(-MaxLag:MaxLag), such that e.g. MIMODisplay(tau,R) can 
%   display the various auto-and cross-correlation functions contained in R.
%
%   Input parameters:
%      X       M x L data matrix
%      MaxLag  maximum lag value calculated
%
%   Output parameter:
%      R       M x M x (2*MaxLag+1) polynomial covariance matrix
%      tau     2*MaxLag+1 index vector for lag values

% S Weiss, Univ of Strathclyde, 20/6/2006

% parameters and initialisation 
[M,L] = size(X);
if MaxLag>=L,
   error('maximum lag cannot exceed length of data');
end;
R = zeros(M,M,2*MaxLag+1);
Lm = L - 2*MaxLag;

% correlation 
Xref = X(:,MaxLag+1:MaxLag+Lm);     % fixed component
for tau = 1:2*MaxLag+1,
  Xshift = X(:,tau:tau+Lm-1);       % delayed component
  R(:,:,2*MaxLag+2-tau) = Xref*Xshift'/Lm; % corr.
end;
tau = (-MaxLag:MaxLag);

% enforce parahermitian property
R = (R + ParaHerm(R))/2;
