function [H2,rho] = PUPolyMatTrim(H1,gamma)
%[H2,rho] = PUPolyMatTrim(H1,gamma)
% 
%  H2=PUPolyMatTrim(H1) trims the time dimension of a paraunitary matrix by
%  removing outer matrix coefficients. The trimming is performed from both
%  ends, and will remove the matrix coefficient with the smallest Frobenius
%  norm, such that the total energy of the trimmed parts does not exceed  
%  one permille of the total energy in H1.
% 
%  H2=PUPolyMatTrim(H1,gamma) trims such that the ratio between the removed 
%  and total energy of R1 is less than gamma, with 0<=gamma<1.
%
%  [H2,rho]=PUPolyMatTrim(H1) or [H2,rho]=PUPolyMatTrim(H1,r)
%  additionally returns the error in paraunitarity --- the norm of 
%  (H2*~H2-I) --- after trimming in rho.
%  
%  This function is based on the 'trim' operation described in [1,2].
%
%  Input parameters:
%     H1      paraunitary matrix
%     gamma   maximum ratio between the maximum energy removed at outer lags 
%             and the total energy 
%             default: 1/1000
%  
%  Output parameters:
%     H2      trimmed parahermitian matrix
%     rho     error in paraunitarity
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
E = gamma*PolyMatNorm(H1);
[M,N,L] = size(H1);

% trim from sequentially from front or back, where ever the matrix with the 
% smallest norm can be found, until the desired amount of energy has been 
% removed
Etrim = 0;
Index_front = 1; Index_back = L;
Efront = norm(H1(:,:,Index_front),'fro').^2;
Eback = norm(H1(:,:,Index_back),'fro').^2;
while (Etrim + min([Efront Eback])) <= E,
   if Eback > Efront,
      % trim at the front
      Etrim = Etrim + Efront;
      Index_front = Index_front + 1;
      Efront = norm(H1(:,:,Index_front),'fro').^2;
   else
      % trim at the back
      Etrim = Etrim + Eback;
      Index_back = Index_back - 1;
      Eback = norm(H1(:,:,Index_back),'fro').^2;
   end;
end;    

% output
H2 = H1(:,:,Index_front:Index_back);

% prepare additional output parameter if error in paraunitarity is requested
if nargout>1,
   HH2 = PolyMatConv(H2,ParaHerm(H2));
   L2 = (size(HH2,3)+1)/2;
   HH2(:,:,L2) = HH2(:,:,L2)-eye(N);
   rho = PolyMatNorm(HH2);
end;
