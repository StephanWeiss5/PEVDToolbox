function C = CodingGain(A);
%C = CodingGain(A);
%
%  CodingGain(vars) returns the coding gain of a subband coding system of M 
%  channels, where vars is an M-element vector of powers. 
%
%  CodingGain(R) returns the coding gain of a subband coding system of M 
%  channels, R is an MxMxL matrix expressing the polynomial covariance matrix
%  of the subchannels. 
%
%  The coding gain is the ratio between arithmetic and geometric mean of the 
%  subband variances.
%
%  Input parameter:
%       vars or R  vector of subband powers / variances, or
%                  polynomial covariance matrix
%
%  Output parameter:
%       C          coding gain

% S. Weiss, 25/8/2014

[M,N,L] = size(A);
if (L==1) && ((M==1)||(N==1)),
   J = max(M,N);
   C = (sum(A)/J) / ((prod(A)).^(1/J));
elseif M==N,
   L2 = (L+1)/2;
   A = diag(A(:,:,L2));
   C = (sum(A)/M) / ((prod(A)).^(1/M));
else
   error('input to CodingGain() does not satisfy the required dimensions');
   C=-1;
end;
 
