function Fnorm = PolyMatNorm(H,SpecString);
%N = PolyMatNorm(H,spec);
%
%  PolyMatNorm(H) returns the generalised square Frobenius norm of all matrix 
%  elements in the polynomial matrix H(z) represented by H. This norm represents 
%  the the sum over all modulus-squared matrix elements, both spatially and 
%  temporally.
% 
%  PolyMatNorm(H,'OnDiag') only returns the sum over of all modulus-squared 
%  diagonal elements. 
%  
%  PolyMatNorm(H,'OffDiag') only returns the sum over all off-diagonal modulus
%  squared elements of H(z).
%
%  PolyMatNorm(H,'Full') is the default and considers all matrix elements.
% 
%  Input parameters:
%     H         MIMO system or polynomial matrix
%     spec      optional parameter:
%               'OnDiag' count on-diagonal elements only
%               'OffDiag' count off-diagonal elements only
%               'Full' count all elements
%               default: 'Full'
%  
% Output parameters:
%     N         norm
  
%  S. Weiss, UoS, 5/10/2005
  
if nargin==1,
   SpecString='Full';
end;
[M,N,L] = size(H);

Fnorm = 0;
if strcmp(SpecString,'Full')==1,
   for l = 1:L,
      A = H(:,:,l);
      Fnorm = Fnorm + sum(sum(A.*conj(A)));
   end;
else
   MNmin = min(M,N);
   Mask = zeros(M,N);
   Mask(1:MNmin,1:MNmin) = eye(MNmin);
   if strcmp(SpecString,'OnDiag')==1,
      Mask = zeros(M,N);
      Mask(1:MNmin,1:MNmin) = eye(MNmin);
   elseif  strcmp(SpecString,'OffDiag')==1,        
      Mask = ones(M,N);
      Mask(1:MNmin,1:MNmin) = Mask(1:MNmin,1:MNmin) - eye(MNmin);
   else
      error('norm option in PolyMatNorm() not defined');
   end;
   for l = 1:L,
      A = H(:,:,l).*Mask;
      Fnorm = Fnorm + sum(sum(A.*conj(A)));
   end;
end;  
