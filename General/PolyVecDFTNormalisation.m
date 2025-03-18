function qf = PolyVecDFTNormalisation(af,Smoothing);
% qf = PolyVecDFTNormalisation(af,Smoothing);
%
% If af is an MxNfft matrix that in its rows represents the DFT of an M-element
% polynomial vector, the return argument qf is the DFT of a vector of analytic 
% functions, approximated to order (Nfft-1), that is a version of af normalised 
% on the unit circle.
% qf = PolyVecDFTNormalisation(af,Smoothing) with Smoothing='on' performs
% additional phase smoothing to the normalised vector. 
%
% The algorithm is based on a draft on Gram-Schmidt orthonormalisation [1], and 
% utilises a bin-shift method detailed in [2]. A smoothing method is utilised
% from [3]; this is a computationally slow and costly procedure, particularly 
% for large DFT sizes!
%
% Input parameter:
%    af         MxNfft  DFT sample points of an M-element polynomial vector
%    Smoothing  'on' for additional phase smoothing (optional)
%
% Output parameter:
%    qf         MxNfft  DFT sample points of M-element normalised vector
%
% References:
% [1]  F.A. Khattak, I.K. Proudler, S.J. Schlecht, and S. Weiss: "Modified Gram-
%      Schmidt Orthonormalisation for Matrices of Analytic Functions," to be
%       submitted.
% [2]  F.A. Khattak, I.K. Proudler, and S. Weiss: "Scalable Analytic Eigenvalue 
%      Extraction Algorithm," submitted to IEEE Access, 2024.
% [3]  S. Weiss, I.K. Proudler, F.K. Coutts, and F.A. Khattak: "Eigenvalue 
%      decomposition of a parahermitian matrix: extraction of analytic 
%      eigenvectors," IEEE Transactions on Signal Processing, vol. 71, pp. 
%      1642-1656, April 2023.

% Stephan Weiss, UoS, 12/8/24

%------------------------------------------------------------------------------
% parameters
%------------------------------------------------------------------------------
[M,K] = size(af);              % spatial dimension and DFT length
eps0 = 1e-7;                   % threshold for "zero norm"
                               % still to be done: should be relative to the largest bin norm
ShiftStep = 0.05;              % shift of any zero bins relative to bin size
QuietMode=0;                   % '0' for intermediate results, '1' for none

%------------------------------------------------------------------------------
% take K-point DFT and move bins where the vector norm is zero
%------------------------------------------------------------------------------
NormAf2 = sum(abs(af).^2);     % squared norm in each bin
ZeroBinIndices = find(NormAf2<eps0);
if QuietMode==0, disp(sprintf('number of bins with zero norm: %d',length(ZeroBinIndices))); end;
% shift frequencies in bins with zeros until norm is finite
ShiftedFreqs = (ZeroBinIndices-1)/K;         % normalised freq. [0...1]
if length(ZeroBinIndices)>0,
  a = ifft(af,K,2);            %   time domain samples needed for interpolation
  ShiftedFreqs = ShiftedFreqs + ShiftStep/K;
  for i = 1:length(ZeroBinIndices),
     Crit = 0; NShifts = 1;
     while Crit==0,
        af(:,ZeroBinIndices(i)) = a*exp(-1i*2*pi*ShiftedFreqs(i)*(0:K-1)');
        NormAf2(ZeroBinIndices(i)) = sum(abs(af(:,ZeroBinIndices(i))).^2);
        if NormAf2(ZeroBinIndices(i))>eps0,      % norm no longer zero
           Crit=1;
        else                                     % keep shifting further
           ShiftedFreqs(i) = ShiftedFreqs(i)+ShiftStep/K;
           NShifts = NShifts+1;
        end;
     end;
     if QuietMode==0, disp(sprintf('bin %d shifted %d times',[ZeroBinIndices(i), NShifts])); end;
     % still need to check that this does not infringe on the next bin
     if NShifts*ShiftStep>=1, error('conflicting bin shift --- either vector norm is too small, or the threshold to detect zeros is too stringent'); end;
  end;
end;

%------------------------------------------------------------------------------
% bin-wise normalisation
%------------------------------------------------------------------------------
NormFactor = 1./sqrt(NormAf2);
Afnorm = af.*NormFactor(ones(1,M),:);

%------------------------------------------------------------------------------
% determine zero crossings 
%------------------------------------------------------------------------------
HermAngle = zeros(1,K);
HermAngle(1) = angle(Afnorm(:,end)'*Afnorm(:,1));
for k = 2:K,
   HermAngle(k) = angle(Afnorm(:,k-1)'*Afnorm(:,k));
end;   
SignCorrection = -mod(cumsum(HermAngle>pi/2),2);      % '-1' means a sign change
dummy = find(SignCorrection==0);
SignCorrection(dummy) = ones(1,length(dummy));        % '1' means no sign change

%------------------------------------------------------------------------------
% fractional delay if number of zero crossings is odd 
%------------------------------------------------------------------------------
if SignCorrection(K)==(-1),
   if QuietMode==0, disp('a sign correction is needed'); end;
   Correction = SignCorrection.*exp(-1i*pi*(0:K-1)/K);
   % fix at shifted frequencies
   Correction(ZeroBinIndices) = SignCorrection(ZeroBinIndices).*exp(-1i*pi*ShiftedFreqs);
else
   Correction = SignCorrection;      
end;
% apply correction
Afnorm = Afnorm.*Correction(ones(1,M),:);   

%------------------------------------------------------------------------------
% interpolate bins with (formerly) zero norm
%------------------------------------------------------------------------------
if length(ZeroBinIndices)>0,                          % skip if there are no zeros
   % matrices for interpolation (see Faizan's IEEE Access submission)
   NonZeroBinIndices = setdiff((1:K),ZeroBinIndices);
   W = dftmtx(K)/sqrt(K);
   Wu = W(NonZeroBinIndices,:);
   Wuperp = W(ZeroBinIndices,:);
   Wn = exp(-1i*2*pi*ShiftedFreqs'*(0:K-1))/sqrt(K);
   A = Wn*Wu'; 
   Binv = inv(Wn*Wuperp');
   % perform interpolation to replace shifted bins
   for m = 1:M,
      Afnorm(m,ZeroBinIndices) = ( [-Binv*A Binv]*Afnorm(m,[NonZeroBinIndices ZeroBinIndices]).' ).';
   end;
end;            
        
%------------------------------------------------------------------------------
% potential smoothing 
%------------------------------------------------------------------------------
P = 3;                          % derivative order applied for smoothing
if nargin>1,
   if strcmp(Smoothing,'on')==1,  
       disp('perform additional phase smoothing');
      qf = PhaseSmoothing(Afnorm,P);
   else
      qf = Afnorm;
   end;
else
   qf = Afnorm;
end;         
