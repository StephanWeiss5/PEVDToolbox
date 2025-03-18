function R = MuxPolyCovMat(r,M);
%R = MuxPolyCovMat(r,M);
%  
%  MuxPolyCovMat(r,M) returns the polynomial covariance matrix arising for a
%  signal characterised by an autocorrelation sequence r that is demulti-
%  plexed into M channels.
%
%  This problem arises e.g. in optimal subband coding [1]. The returned 
%  matrix R is pseudo-circulant.
%
%  Input parameters:
%      r       autocorrelation sequence of input process
%      M       number of subchannels for demultiplexing
%
%  Output parameters:
%      R       MxMxL polynomial covariance matrix
%
%  Reference:
%
%  [1] S. Redif, J.G. McWhirter, and S. Weiss, "Design of FIR Paraunitary 
%      Filter Banks for Subband Coding Using a Polynomial Eigenvalue Decompo- 
%      sition," IEEE Transactions on Signal Processing, vol. 59, no. 11, 
%      pp. 5253-5264, Nov 2011.

% S Weiss, UoS, 16/9/2004
% updated, S Weiss, Strathclyde, 27/8/14

[N1,N2] = size(r);
if N1==1,
   r=r.';
   L=N2;
else
   L=N1;
end;
if mod(L,2)~=1,
   error('MuxPolyCovMat() fed with an ACS of even length');
end;  
Lr = (L+1)/2;             % number of lags in r to one side
r  = r(Lr:end);
LR = ceil((Lr+M-1)/M);    % number of lags in R to one side

Rfull = toeplitz([r; zeros(2*M,1)]);
R(:,:,LR) = Rfull(1:M,1:M);
for i = 1:(LR-1),
   R(:,:,LR+i) = Rfull(1:M,i*M+1:(i+1)*M);
   R(:,:,LR-i) = R(:,:,LR+i)';
end;
  
