function [Qstar,ePU,Nfft] = NearPUProcrustes(A,NDFTmax,PUerror);
% NearPUProcrustes(A,NFTmax,PUerror); 
%
% Qstar = NearPUProcrustes(A,NFTmax,PUerror) takes a near-paraunitary matrix A
% and calculates a bin-wise Procrustes solution in Qstar. The DFT size is iter-
% atively increased until either a threshold for the error in paraunitary, 
% PUerror, is met, or a maximum DFT size NDFTmax is reached.
% 
% Near PU means that A(z) should only possess real positive singular values on the
% unit circle.
%
% [Qstar,ePU,NDFT] = NearPUProcrustes(A,NFTmax,PUerror) additionally returns the
% reached error in paraunitarity as the required DFT length.
%
% Input parameters:
%    A          MxMxL matrix to be approximated (polynomial coefficients)
%    NDFTmax    maximum DFT length where iterations are stops
%    PUerror    threshold for error in paraunitarity
%
% Output parameters:
%    Qstar      MxMxL2 Procrustes solution (Polynomial coefficients)
%    ePU        paraunitarity error achieved by the solution
%    Nfft       DFT length at which iterations terminated

% S. Weiss, UoS, 1/2/2025

%----------------------------------------------------------------------------
%  parameters and initialisations
%----------------------------------------------------------------------------
[M,~,L] = size(A);
Nfft = 2^(ceil(log2(L)));

Af = fft(A,Nfft,3);
Qf = zeros(M,M,Nfft);
for k = 1:Nfft,
   [u,~,v] = svd(squeeze(Af(:,:,k)));
   Qf(:,:,k) = u*v';
end;
         
%----------------------------------------------------------------------------
%  iteration
%----------------------------------------------------------------------------
IterCrit = 1; 
while IterCrit == 1,
   Nfft = Nfft*2;
   % calculate bins that had previously not been considered
   Amod = A;
   for ll = 1:L, Amod(:,:,ll) = Amod(:,:,ll)*exp(-1i*2*pi/Nfft*(ll-1)); end;
   Afrem = fft(Amod,Nfft/2,3);       % contains the remaining, unconsidered DFT bins
   % only calculate bin-wise SVDs in bins that had not been evaluated yet
   Qfnew = zeros(M,M,Nfft);
   Qfnew(:,:,1:2:end) = Qf;
   for k = 1:1:Nfft/2,
     [u,~,v] = svd(Afrem(:,:,k));
     Qfnew(:,:,2*k)=u*v';
   end;
   
   % perform time-domain reconstruction
   Qnew = circshift(ifft(Qfnew,Nfft,3),[0 0 Nfft/2]);   
   % trim    
   p = squeeze(sum(sum(abs(Qnew).^2,1),2));                     % energy per coefficient
%   plot(10*log10(p)); drawnow; pause;
   NonZeroIndicesPre = find(10*log10(p)>-300);                  % find range of non-zero coeff.
   NStart = NonZeroIndicesPre(1); NEnd = NonZeroIndicesPre(end);
   Qstar = Qnew(:,:,NStart:NEnd);
   % check error in paraunitarity
   ePU = PUMismatch(Qstar);
    
   if (ePU<PUerror)||(Nfft>NDFTmax/2),
      IterCrit=0;
   else
      Qf = Qfnew;
   end;   
end;      
   
   
   
   
   
   
