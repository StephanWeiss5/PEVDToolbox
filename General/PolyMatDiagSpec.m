function P = PolyMatDiagSpec(H,Ndft);
%P = PolyMatDiagSpec(R,Ndft);
% 
%   PolyMatDiagSpec(R) calculates the spectra on the main diagonal of a poly-
%   nomial matrix R of dimension MxNxL, evaluated over L DFT bins. The 
%   spectra are returned as the columns of an LxK matrix P, with K=min(M,N).
%
%   PolyMatDiagSpec(R,Ndft) calculates the spectra using Ndft number of bins,
%   whereby Ndft must be greater or equal L. The spectral are returned as 
%   the columns of an (Ndft)xK matrix P.
%
%   Input parameters:
%      R       MxNxL polynomial matrix
%      Ndft    number of DFT bins (options)
%              default: L  
%
%   Output parameter:
%      P       Ndft x min(M,N) matrix of spectra

% S. Weiss, 20/10/2005
  
[M,N,L] = size(H);  
  
if nargin==1,
   Ndft = L;
end;
if Ndft<L,
   warning('parameter Ndft in function MIMODiagSpec() is too small');
end;

MN = min(M,N);
P = zeros(Ndft,MN);
for i = 1:MN,
   P(:,i) = fft(shiftdim(H(i,i,:),2),Ndft);
end;

