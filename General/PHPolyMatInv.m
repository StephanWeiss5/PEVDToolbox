function S = PHPolyMatInv(R,T,MaxIter,theta,epsilon);
%S = PHPolyMatInv(R,T,MaxIter,theta,epsilon);
%
%  PHPolyMatInv() calculates the inverse of parahermitian polynomial matrix as
%  described in [1].
%
%  S = PHPolyMatInv(R) returns the inverse of R. The input R  represents an 
%  MxMx(2L+1) parahermitian matrix R(z) of the form
%     R(z) = RL'z^L + ... + R1' z + R0 + R1 z^{-1} + ... + RLz^{-L}
%  whereby
%     R(:,:,1) = RL';
%     ...
%     R(:,:,L) = R1';
%     R(:,:,L+1) = R0;
%     R(:,:,L+2) = R1;
%     ...
%     R(:,:,2*L+1) = RL;
%  The function returns a paraunitary matrix S(z) such that approximately
%  S(z) R(z) = R(z) S(z) = I. The polynomial matrix S(z) is returned in the 
%  same format at R, with the lag parameter forming a 3rd dimension.
%
%  The inversion is based in a polynomial EVD of R, and an inversion of the
%  polynomial eigenvalues. The PEVD uses the sequential best rotation (SMD, [2])
%  algorithm, and S(z) is determined from a minimum mean square error inverse
%  detailed in [1]. 
%
%  S = PHPolyMatInv(R,T) sets the order of the inverse polynomial eigenvalues 
%  to 2T (default: T=10*L). 
% 
%  The SMD algorithm can be influenced by a number of optional inputs:
%
%  S = PHPolyMatInv(R,T,MaxIter) stops the SMD algorithm after MaxIter 
%  iteration steps (default: MaxIter=100).
%
%  S = PHPolyMatInv(R,T,MaxIter,theta) stops the SMD algorithm either 
%  after MaxIter iteration steps, or once the maximum off-diagonal element 
%  falls below a threshold theta (default: theta=1e-5).
% 
%  S = PHPolyMatInv(R,T,MaxIter,theta,epsilon) curtails the growth of
%  both the parahermitian and paraunitary matrices that arise from the PEVD
%  of R, by removing a proportion epsilon of the energy from these matrices.
%
%  Input parameters
%     R         parahermitian matrix
%     T         length of SISO inverses
%     MaxIter   max. iterations for PEVD
%               (default = 200)
%     theta     threshold for off-diagonal values
%               (default = 0.00001)
%     epsilon   proportion of energy being trimmed to curtail polynomial order
%
%  Output parameter:
%     S         inverse parahermitian matrix
%
%  References:
%  [1] S. Weiss, A. Millar, and R.W. Stewart: "Inversion of Parahermitian 
%      Matrices,"European Signal Processing Conference, Aalborg, Denmark, 
%      pp. 447-451, August 2010.
%  [2] S. Redif, S. Weiss, and J.G. McWhirter, "Sequential Matrix Diagonali-
%      sation Algorithms for Polynomial EVD of Parahermitian Matrices," IEEE
%      Transactions on Signal Processing, 63(1):81-89, January 2015. 

% S.Weiss, UoS, updated 6/6/2015

%------------------------------------------------
%   check for optional parameters
%------------------------------------------------
if nargin<5, epsilon = 0.0; end;
if nargin<4, theta = 0.00001; end;
if nargin<3, MaxIter = 200; end;
if nargin<2, T = 10*size(R,3); end;

%------------------------------------------------
%   polynomial EVD
%------------------------------------------------
[Q,Gamma] = SMD(R,MaxIter,theta,0.0);
GammaInv = DiagInverse(Gamma,T,1e-5);
S = PolyMatConv(ParaHerm(Q),PolyMatConv(GammaInv,Q));



%------------------------------------------------
%  subroutine DiagInverse()
%------------------------------------------------
function Hinv = DiagInverse(D,T,gamma);
% Inversion of a diagonal polynomial matrix D. D does not have
% to be diagonal, but for the calculation of Dinv, only the on-diagonal
% elements of D are considered.
% 
% The inverse is calculated on the index interval [-T;+T]. 
%  
% The diagonal elements of D are truncated if there are tails with
% elements smaller than the variable gamma.
%
% Input parameters:
%    D         KxK polynomial matrix
%    T         length of inverse (or lag of its autocorrelation sequence)
%    gamma     element size for truncation
%
% Output parameters:
%    Dinv      KxK diagonal matrix containing the inverse
%              polynomials of the on-diagonal elements of H.
%  
% S. Weiss, 28/12/2006

[M,N,L] = size(D);
Dinv = zeros(M,M,T);
for m = 1:M,
   d = shiftdim(D(m,m,:));
   ZeroPosition = find(abs(d)>gamma);        % identify zeros
   d = d(ZeroPosition(1):ZeroPosition(end));   % truncate
   hinv = AcsMmseInv(d,T);
   Hinv(m,m,:) = hinv;
end;


%------------------------------------------------
%  subroutine AcsMmseInv()
%------------------------------------------------
function s = AcsMmseInv(r,T);
% Inverts the response r which is assumed to have the
% properties of an autocorrelation sequence. The inverse
% is constrained to the same symmetry conditions as r,
% with a length of 2T+1, i.e. corresponding to an inverse
% autocorrelation sequence with maximum lag +/-T.
%
% Input parameters
%    r     autocorrelation sequence
%    T     maximum lag of inverse
%
% Output parameter
%    s     inverse of r
%
% S. Weiss, univ of Strathclyde, 24/1/2010

%-----------------------------------------------------
% setup matrix equation
%-----------------------------------------------------
T2 = (T+1)/2;
R = toeplitz([r; zeros(T-1,1)],[r(1) zeros(1,T-1)]);
R1 = R(:,1:T2) + fliplr(R(:,T2:end));
R2 = R(:,1:T2) - fliplr(R(:,T2:end));
R_real = [real(R1)  -imag(R2);
        imag(R1)   real(R2)];
d_real = zeros(2*(length(r)+T-1),1);
d_real( (length(r)+T)/2 ) = 1;

%-----------------------------------------------------
% pseudo-inverse
%-----------------------------------------------------
s_real = pinv(R_real)*d_real;

%-----------------------------------------------------
% assemble inverse
%-----------------------------------------------------
s = [s_real(1:T2); zeros(T2-1,1)] + sqrt(-1)*[s_real(T2+1:2*T2); ...
                    zeros(T2-1,1)];
s = s + conj(flipud(s));

