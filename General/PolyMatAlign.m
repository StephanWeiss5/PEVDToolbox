function [Shift,error,A2,B2] = PolyMatAlign(A,B);
%[Shift,error] = PolyMatAlign(A,B);
%
%   PolyMatAlign(A,B) takes input arguments A and B to represent polynomial 
%   matrix A(z) and B(z) of the same dimension but potentially different order.
%   It determines the best delay such that A(z) - B(z)z^{-Shift} is minimised
%   in the least squares sense.
%   
%   [Shift,error] = PolyMatAlign(A,B) returns the determined shift and norm of
%   of the mismatch between the aligned matrices A(z) and B(z)z^{-Shift}.
%
%   [Shift,error,A2,B2] = PolyMatAlign(A,B) additionally returns the aligned 
%   matrices of identical orders in A2 and B2.
%
%   Input parameters:
%      A           MxNxL1 matrix
%      B           MxNxL2 matrix
%   
%   Output parameters:
%      Shift       determined delay between A and B
%      error       least squares mismatch between the aligned matrices
%      A2          MxBxL3 matrix A after alignment
%      B2          MxBxL3 matrix B after alignment

% S. Weiss, 19/2/2023

%----- check arguments
[M1,N1,L1] = size(A); [M,N,L2] = size(B);
if (M~=M1)||(N~=N1), error('dimension mismatch between input arguments'); end;

%------------------------------------------------------------------------------
%   determine shift
%------------------------------------------------------------------------------
%----- elementwise auto-correlation, add moduli
r_ab = zeros(L1+L2-1,1);
for m = 1:M, 
   for n = 1:N,
       r_ab = r_ab + abs(conv(squeeze(A(m,n,:)),flipud(conj(squeeze(B(m,n,:))))));
   end;    
end;

%----- determine delay
TimeScale=(-L2+1:L1-1);
[~,TIndex] = max(r_ab);
Shift = -TimeScale(TIndex);

%------------------------------------------------------------------------------
%   determine mismatch
%------------------------------------------------------------------------------
% compensate for shift
if Shift>0,    % a positive shift means that A needs zero-padding at the front
   dummy = zeros(M,N,L1+Shift);
   dummy(:,:,Shift+1:end) = A;
   A = dummy;
else           % a negative shift means that B needs zeropadding at the front
   dummy = zeros(M,N,L2-Shift);
   dummy(:,:,1-Shift:end) = B;
   B = dummy;
end;    
% compensate for difference in orders
L1 = size(A,3); L2 = size(B,3);
if L1>L2,
   B2 = zeros(M,N,L1);
   B2(:,:,1:L2) = B;
   A2 = A;
else
   A2 = zeros(M,N,L2);
   A2(:,:,1:L1) = A;
   B2 = B;
end;
% mismatch as sum of squared Frobenius norms
error = PolyMatNorm(A2-B2); 
