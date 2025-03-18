function u = PolyVecNormalisation(v,rho,N,epsilon);
%u = PolyVecNormalisation(v,rho); 
%
%  This function normalises a polynomial vector v to unit length, such that
%  u^P(z).u(z) = 1. This is equivalent to normalising a constant vector to
%  a Euclidean length of unity in the non-polynomial case.
%
%  The result vector is trimmed to remove a fraction rho of the energy before
%  returning u, in order to limit the order of this vector.
%
%  Implemented in the DFT domain, the chosen DFT size exceeds 2^N times the 
%  order of v, with a default of N=4. This order increase is for internal 
%  precision, and the result is trimmed subsequently.
%  
%  Input parameter:
%       v           Mx1xL matrix containing a polynomial vector of length L
%       rho         fraction of energy to be trimmed for return variable
%       N           internal order increase is 2^N (default: N=4)
%       epsilon     small constant to avoid division by zero (default 2.2e-16)
%
%  Output parameter:
%       u           Mx1xJ matrix of normalised polynomial vector
%        
        
%  S. Weiss, UoS, 27/7/2022

if nargin<3,
   N = 4;
   epsilon = eps;
elseif nargin<4,
   epsilon=eps;
else
   % do nothing   
end;   
%epsilon = 1e-16;
%epsilon

%------------------------------------------------------------------------------
% determine normalisation function;
%------------------------------------------------------------------------------
rvv = squeeze(PolyMatConv(ParaHerm(v),v));
Lrvv = length(rvv);
Nfft = 2^(ceil(log2(Lrvv))+4);
Rvv = sqrt(abs(fft(rvv,Nfft))); 
%------------------------------------------------------------------------------
% perform normalisation in the DFT domain
%------------------------------------------------------------------------------
U = fft(v,Nfft,3);
for k = 1:Nfft,
   if Rvv(k) <= epsilon
      U(:,:,k) = 0;
   else   
      U(:,:,k) = U(:,:,k)/Rvv(k);
   end;   
end;
u = ifft(U,Nfft,3);

%------------------------------------------------------------------------------
% allow for non-causal solution
%------------------------------------------------------------------------------
dummy = u;
u(:,:,1:Nfft/2) = dummy(:,:,Nfft/2+1:Nfft);
u(:,:,Nfft/2+1:Nfft) = dummy(:,:,1:Nfft/2);
   
%------------------------------------------------------------------------------
% trimming
%------------------------------------------------------------------------------
E = zeros(Nfft,1);
for k = 1:Nfft,
   E(k) = norm(u(:,:,k)).^2;
end;         
Etotal = sum(E);    
% now take iteratively from the front or the end of u until a fraction rho of 
% the energy has been removed 
E = E/Etotal;
Trim=1; StartIndex=1; EndIndex=Nfft; ClippedEnergyFraction=0;
while (Trim==1),
  if E(StartIndex)<E(EndIndex),
     if (ClippedEnergyFraction+E(StartIndex))<=rho,
        ClippedEnergyFraction = ClippedEnergyFraction+E(StartIndex);
        StartIndex=StartIndex+1;
     else
        Trim=0;
     end;
  else
     if (ClippedEnergyFraction+E(EndIndex))<=rho,
        ClippedEnergyFraction = ClippedEnergyFraction+E(EndIndex);
        EndIndex=EndIndex-1;
     else
        Trim=0;
     end;
  end;
end;
u = u(:,:,StartIndex:EndIndex);                


