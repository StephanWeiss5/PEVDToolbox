function [H,D,F] = GenerateMIMOSources(L,P,M,K,gamma,Mode)
%[H,D,F] = GenerateMIMOSources(L,P,M,K,gamma,Mode)
%
%  GenerateMIMOSources(L,P,M,K) produces a source model for M-array data 
%  emanating from L independent Gaussian sources spectrally shaped by Pth 
%  order moving average (MA, i.e. finite impulse response) innovation 
%  filters. The source signals are then mixed by a random MxL paraunitary 
%  matrix of order K.
%
%  GenerateMIMOSources(L,P,M,K,gamma) operates the same but the additional 
%  input parameter gamma limits the radii of zeros in the MA innovation 
%  filters. This restricts the dynamic range of the generated spectra.
%
%  A flow graph of this source model is given in [1]. In the context of 
%  multichannel coding, M=L is required to avoid a division by zero, as the 
%  geometric mean of variances is zero for M>L. The purpose of limiting the 
%  dynamic range of source PSDs is driven by the same aim, and keeps the 
%  problem numerically tractable.
%
%  The source model does currently not permit the overdetermined case where 
%  the number of sources exceeds the number of sensors.
%
%  [H,D,F]=GenerateMIMOSources(L,P,M,K,gamma) returns the paraunitary 
%  mixing matrix H(z) of dimension MxL and order K in H; D represents a 
%  diagonal parahermitian matrix D(z) containing the power spectra of the L 
%  sources prior to convolutive mixing. Minimum phase filters of order P 
%  that generate such PSDs from unit variance complex Gaussian white noise 
%  are contained in the columns of the matrix F.
%
%  [H,D,F]=GenerateMIMOSources(L,P,M,K,gamma,Mode) returns a real-valued
%  source model for Mode='real', and a complex valued one by default.
%
%  Input parameters
%     L       number of independent sources
%     P       order of MA innovation filter
%     M       dimension of array generated (M>=L)
%     K       order of paraunitary mixing matrix
%     gamma   max radius for zeros in the innovation filters
%           (default value is 1)
%     Mode    'real' for readl valued, 'complex valued' (default)
%
%  Output parameters
%     H       paraunitary matrix
%     D       diagonal, parahermitian matrix
%     F       columns are the innovation filters, whose 
%             ACFs generate the diagonal entries of D.
%
%  Reference:
%
%  [1] S. Redif, S. Weiss and J.G. McWhirter, "Sequential Matrix 
%      Diagonalisation Algorithms for Polynomial EVD of Parahermitian 
%      Matrices," IEEE Transactions on Signal Processing, 63(1):81-89, 
%      January 2015.

% Stephan Weiss, August 2012
% real valued option, S. Weiss, Oct 2018

%--------------------------------------------------------
%   Input Parameter Check
%--------------------------------------------------------
if nargin==4,
   gamma=1;
end;
if M < L,
  error('more sources than array elements');
end;
RealMode = 0;
if nargin==6,
  if Mode=='real',
     RealMode=1;
  end;
end;  

%--------------------------------------------------------
%   Paraunitary Matrix
%--------------------------------------------------------
H = eye(M);  
for k = 1:K,
  if RealMode==1, 
     v = randn(M,1);
  else
     v = randn(M,1) + sqrt(-1)*randn(M,1);
  end;
  v = v./norm(v,2);
  Hnew = zeros(M,M,2);
  if RealMode==1,
    [Q,~] = qr(randn(M,M));
  else  
    [Q,~] = qr(randn(M,M) + sqrt(-1)*randn(M,M));  
  end;  
  Hnew(:,:,1) = v*(v'*Q);
  Hnew(:,:,2) = Q - Hnew(:,:,1);
  H = PolyMatConv(H,Hnew);
end;

%--------------------------------------------------------
%   Source Innovation Filters
%--------------------------------------------------------
if RealMode==1,
  P2 = floor(P/2);
  P = 2*P2;
  Fz = (rand(P2-1,L)*gamma).*exp(sqrt(-1)*2*pi*rand(P2-1,L));  % zeros
  Fz = [Fz; conj(Fz)];
else    
  Fz = (rand(P-1,L)*gamma).*exp(sqrt(-1)*2*pi*rand(P-1,L));  % zeros
end;
F = zeros(P,L);
for l = 1:L,
   F(:,l) = poly( Fz(:,l)).';
end;
F(:,1) = F(:,1)./norm(F(:,1),2);
Ffd = abs(fft(F,1024,1));
for l = 1:L-1,
   alpha = min(Ffd(:,l)./Ffd(:,l+1));
   Ffd(:,l+1) = alpha*Ffd(:,l+1);
   F(:,l+1) = alpha*F(:,l+1);
end;

%--------------------------------------------------------
%   Output Parameters
%--------------------------------------------------------
D1 = zeros(L,L,P);
for l = 1:L,
  for p = 1:P,
    D1(l,l,p) = F(p,l);
  end;
end;
D = PolyMatConv(D1,ParaHerm(D1));
H = H(:,1:L,:);



