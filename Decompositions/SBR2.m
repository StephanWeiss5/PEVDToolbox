function [H,Gamma] = SBR2(R,maxiter,epsilon,Mu,vers);
%[H,Gamma] = SBR2(R,maxiter,epsilon,mu,vers);
%
%  Polynomial matrix eigenvalue decomposition (PEVD) algorithm factorising
%  a parahermitian matrix represented by R using the second order sequential
%  best rotatio (SBR2) algorithm [1].
%
%  [H,Gamma]=SBR2(R) takes as input an MxMx(2L+1) matrix R representing a 
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
%  SBR2 will strive of diagonalise and spectrally majorise Gamma(z). The 
%  algorithm stops either after a fixed number of iterations or once a 
%  threshold for the maximum absolute value of off-diagonal elements has 
%  been reached. Default values are outlined below.
%
%  [H,Gamma] = SBR2(R,maxiter) stops after maxiter iterations. The default
%  value for the optional parameter maxiter is 400.
%
%  [H,Gamma] = SBR2(R,maxiter,epsilon) stops either after maxiter iterations
%  or once the maximum absolute off-diagonal element has a value smaller than
%  epsilon.
%
%  [H,Gamma] = SBR2(R,maxiter,epsilon,mu) additionally performs a truncation
%  of the parahermitian matrix at every iteration such that matrices at outer
%  lags containing a mu-th of the total power are truncated. This stops
%  unnecessary growth of the parahermitian matrix, and therefore keeps 
%  computational complexity down and Gamma of sufficiently low order.
%
%  [H,Gamma] = SBR2(R,maxiter,epsilon,mu,vers) take the optional input vers
%  to switch between SBR2 (vers='SBR2') and the coding-gain optimised SBRC2 
%  (vers='SBR2C').
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
%     vers      SBR2 version ('SBR2' or 'SBR2C')
%               default: 'SBR2C'
%
%  Output parameters:
%     H         paraunitary matrix
%     Gamma     approximately diagonalised and spectrally majorised para-
%               hermitian matrix
%
%  Please not that the SBR2 algorithm is subject of a patent owned by 
%  QinetiQ. QinetiQ has given permission for its free use in university 
%  research, but at present its use in a commercial application without 
%  license from QinetiQ would be an infringement of the patent.
%
%  References:
%
%  [1] J.G. McWhirter, P.D. Baxter, T. Cooper, S. Redif, and J. Foster, "An 
%      EVD Algorithm for Para-Hermitian Polynomial Matrices," IEEE Trans-
%      actions on Signal Processing, vol. 55, no. 5, pp. 2158-2169, May 2007.
%
%  [2] S. Redif, J.G. McWhirter, and S. Weiss, "Design of FIR Paraunitary 
%      Filter Banks for Subband Coding Using a Polynomial Eigenvalue Decompo- 
%      sition," IEEE Transactions on Signal Processing, vol. 59, no. 11, 
%      pp. 5253-5264, Nov 2011.
%
%  Acknowledgement to the appropriate paper should be given if this function 
%  is utilised for academic output. 
