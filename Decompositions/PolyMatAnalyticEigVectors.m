function [Q,chi] = PolyMatAnalyticEigVectors(R,Lambda,Nmax,Thresh1,Thresh2);
%[Q,chi] = PolyMatAnalyticEigVectors(R,Lambda,Nmax,Th1,Th2)
%
%   Extract approximate analytic eigenvectors from a parahermitian matrix R 
%   using the approach described in [4], returning analytic eigenvalues in 
%   Lambda. These are used to order the eigenvectors in each bin, extract 1D 
%   eigenspaces across algebraic multiplicitities, and obtain an analytic 
%   eigenvector in each such 1d subspace through finding the phase in each bin
%   that creates the smoothes possible function, or shortest time-domain 
%   support.
%
%   Input parameters
%       R        MxMxN space-time covariance matrix (time-domain)
%       Lambda   MxL analytic eigenvalues (time domain, one per row)
%       Nmax     maximum FFT length for iteration (default 256)
%       Th1      threshold for orthonormality (default 5e-5)
%       Th2      threshold for trimming (default 1e-4)
%
%   Output parameters
%       Q        MxMx? matrix of eigenvectors (time domain)
%       chi      error in paraunitarity
%
%  References
%  [1] S Weiss, J Pestana and IK Proulder: "On the existence and uniqueness of
%      the eigenvalue decomposition of a parahermitian matrix," IEEE Trans. on
%      Signal Processing, 66(10):2659-2672, May 2018.
%  [2] S Weiss, J Pestana, IK Proulder, and FK Coutts: "Corrections to `On the 
%      existence and uniqueness of the eigenvalue decomposition of a para-
%      hermitian matrix'," IEEE Trans. on Signal Processing, 66(23):6325-6327, 
%      Dec. 2018.
%  [3] S Weiss, IK Proulder, and FK Coutts: "Parahermitian matrix eigenvalue 
%      decomposition: extraction of analytic eigenvalues," IEEE Trans. on 
%      Signal Processing, 69:722-737, Jan. 2021.
%  [4] S Weiss, IK Proulder, FK Coutts, and F Khattak: "Parahermitian matrix eigenvalue 
%      decomposition: extraction of analytic eigenvectors," IEEE Trans. on 
%      Signal Processing, submitted, Feb. 2022.
