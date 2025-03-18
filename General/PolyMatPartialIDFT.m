function At = PolyMatPartialIDFT(Af,Nfft,Bins);
% At = PolyMatPartialIDFT(Af,Nfft,Bins);
%
% PolyMatPartialIDFT() performs an IDFT from only a subset of frequency bins. 
% For the Fourier representation of a polynomial matrix Af, the index vector 
% Bins describes which of a total of Nfft bins are present in Af. For the
% reconstruction to be exact, the orginal FFT has to be at least twice as long
% as the temporal support of Af.
%
% Input parameters:
%     Af               MxMxL matrix in the Fourier domain
%     Nfft             FFT size (must be even)
%     Bins             index into the bins that are present in Af
%
% Output parameter:
%     At               MxMx(Nfft/2)  time domain reconstruction

% S. Weiss, 27/9/23, based on code by Sebastian Schlecht

T = dftmtx(Nfft);
M = size(Af,1);
Tinv = pinv(T(Bins,1:Nfft/2));
At = zeros(M,M,Nfft/2);
for m = 1:M,
   for n = 1:M,
      dummy = Tinv*squeeze(Af(m,n,:));
      At(m,n,:) = dummy; 
   end;
end;      
