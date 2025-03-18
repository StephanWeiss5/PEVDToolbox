function [L_analytic,L_permutation,EVPrecision,xi2] = PolyMatAnalyticEigValues(R,Nmax,PrecLimit);
%PolyMatAnalyticEigValues(R);
%  [L,Perm,prec,xi2] = PolyMatAnalyticEigValues(R) returns the extracted analytic
%  eigenvalues L of the polynomial matrix R [1,2] as described in [3]. The algorithm
%  operated in discrete frequency bins and iterates until a defined bound for xi2 
%  (see [3]) is reached, or the FFT length exceeds Nmax. Compared to a spectrally 
%  majorised ordering, the permutations of the eigenvalues in each frequency bin 
%  are returned in Perm. An approximate approximation error is returned in EVPrecision,
%  which is a mix of truncation and time-domain aliasing. 
%
%  Input parameters
%     R         MxMxL parahermitian matrix
%     Nmax      maximum FFT length (optional; default 2^10)
%     PrecLimit precision limit for iteration
%
%  Output parameters
%     L         MxK matrix containing the extracted analytic eigenvalues in
%               its rows
%     Perm      Permutations compared to the spectrally majorised solution
%     prec      time domain precision
%     xi2       convergence/divergence metric
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
%      Signal Processing, 69:722-737, Jan. 2021
