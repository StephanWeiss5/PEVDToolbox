function [H,Gamma] = SMD(R,maxiter,epsilon,Mu,SFlag);
%[H,Gamma] = SMD(R,maxiter,epsilon,Mu,vers);
%
%  Polynomial matrix eigenvalue decomposition (PEVD) algorithm factorising
%  a parahermitian matrix represented by R using the sequential matrix
%  diagonalisation (SMD) method [1].
%
%  [H,Gamma]=SMD(R) takes as input an MxMx(2L+1) matrix R representing a 
%  parahermitian matrix R(z) of the form
%     R(z) = RL'z^L + ... + R1' z + R0 + R1 z^{-1} + ... + RLz^{-L}
%  whereby
%     R(:,:,1) = RL';
%     ...
%     R(:,:,L) = R1';
%     R(:,:,L+1) = R0;
%     R(:,:,L+2) = R1;
%     ...
%     R(:,:,2*L+1) = RL;
%  The function returns a paraunitary matrix H(z) in H which creates an 
%  approximately diagonalised parahermitian Gamma(z)
%     Gamma(z) = H(z) R(z) H~(z).
%  The format of Gamma representing Gamma(z) is analogously to R above. For
%     H(z) = H0 + H1 z^{-1} + H2 z^{-2} + ...
%  the returned parameter H is  
%     H(:,:,1) = H0;
%     H(:,:,2) = H1;
%     H(:,:,3) = H2;
%      ...
%  SMD will strive of diagonalise and spectrally majorise Gamma(z). The 
%  algorithm stops either after a fixed number of iterations or once a 
%  threshold for the maximum absolute value of off-diagonal elements has 
%  been reached. Default values are outlined below.
%
%  [H,Gamma] = SMD(R,maxiter) stops after maxiter iterations. The default
%  value for the optional parameter maxiter is 400.
%
%  [H,Gamma] = SMD(R,maxiter,epsilon) stops either after maxiter iterations
%  or once the maximum absolute off-diagonal element has a value smaller than
%  epsilon.
%
%  [H,Gamma] = SMD(R,maxiter,epsilon,mu) additionally performs a truncation
%  of the parahermitian matrix at every iteration such that matrices at outer
%  lags containing a mu-th of the total power are truncated. This stops
%  unnecessary growth of the parahermitian matrix, and therefore keeps 
%  computational complexity down and Gamma of sufficiently low order.
%
%  [H,Gamma] = SMD(R,maxiter,epsilon,mu,vers) take the optional input vers
%  to switch between SMD (vers='SMD') and the coding-gain optimised maximum 
%  element SMD (MESMD, vers='MESMD').
%
%  Input parameters:
%     R         polynomial covariance matrix
%     maxiter   maximum number of iterations (optional)
%               default: 400  
%     epsilon   stop if largest absolute off-diag element is smaller than 
%               epsilon (optional) 
%               default: 0.0001  
%     mu        power ratio in tail of polynomial matrix to be truncated at 
%               every iteration (optional)
%               default: 0.0 (only truncating true zeroes)
%     vers      SMD version ('SMD' or 'MESMD')
%               default: 'SMD'
%
%  Output parameters:
%     H         paraunitary matrix
%    Gamma      polynomial covariance matrix
%
%  Reference:
%
%  [1] S. Redif, S. Weiss, and J.G. McWhirter, "Sequential Matrix Diagonali-
%      sation Algorithms for Polynomial EVD of Parahermitian Matrices," IEEE
%      Transactions on Signal Processing, 63(1):81-89, January 2015. 
%
%  Please acknowledge this paper if this function is utilised for academic output. 
