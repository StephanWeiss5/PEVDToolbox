function R0 = ParaHermFourierCoeff(R,f0);
% R0 = ParaHermFourierCoeff(R,f0);
%
%  Given a MxMxL parahermitian matrix in the time domain, this function 
%  evaluates the Fourier coefficient matrix at normalised frequency f0.
%  The normalisation is such that f0 = 0 is DC, and f0 = 1 the sampling 
%  rate.
%
%  Input parameters:
%      R      MxMxL parahermitian matrix in the time/coefficient domain
%      f0     normalised frequenc [0...1]
%
%  Output parameter:
%      R0     MxM Fourier coefficient matrix at norm. freq. f0

% S. Weiss, 18/4/2021

% parameters
[M,~,L] = size(R);
Maxtau = (L-1)/2;

% Fourier coefficient
R0 = zeros(M,M);
for m = 1:M,
   for mu = 1:(m-1),
      R0(mu,m) = squeeze(R(mu,m,:)).'*exp(-sqrt(-1)*2*pi*f0*(-Maxtau:Maxtau)');
   end;
   R0(m,m) = .5*squeeze(R(m,m,:)).'*exp(-sqrt(-1)*2*pi*f0*(-Maxtau:Maxtau)');
end;
R0 = R0+R0';

