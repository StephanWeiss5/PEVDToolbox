function Y = MIMOFilter(H,X);
%Y = MIMOfilter(H,X);
%
%  Performs a filtering of an N-channel input with 
%  a KxN channel MIMO system.
%  
%  Input parameters:
%       H     K x N x L MIMO system matrix
%                 K   output dimension
%                 N   input dimension
%                 L1   lenght of FIR filters
%       X     N x L2   N channel input of length L2
%
%  Output parameter:
%       Y     K x L2   K channel output of length L2
  
%  S Weiss, Univ of Southampton, 15/7/2004
  
[M,N1,L1] = size(H);
[N2,L2] = size(X);
if N2 ~= N1,
  error('MIMOfilter: dimensions do not agree');
end;
Y = zeros(M,L2);
for m = 1:M,
  for n = 1:N1,
    Y(m,:) = Y(m,:) + filter(shiftdim(H(m,n,:)),1,X(n,:));  
  end;  
end;

