function H = PUPolyMatRand(M,L,N,mode);
%H = PUPolyMatRand(M,L,N,mode);
%
%  H = PUPolyMatRand(M,L) generates a random MxM paraunitary matrix of order L.
%  The matrix is assembled from L first order elementary paraunitary matrices
%  described in [1],
%
%  V(z) = I - u u' + u u' z^{-1}
% 
%  where u is an arbitrary unit norm vector. The resulting paraunitary matrix 
%  H(z) is the product of L such elementary paraunitary matrices, s.t.
%    H(z) = H0 + H1 z^{-1} + H2 z^{-2} + ... + HL z^{-L}
%  is represented in a 3-dimensional matrix as
%    H(:,:,1) = H0;
%    H(:,:,2) = H1;
%    H(:,:,3) = H2;
%      ...
%    H(:,:,L) = HL;
%  Paraunitarity means that both PolyMatConv(H,ParaHerm(H)) and 
%  PolyMatConv(ParaHerm(H),H) will result in an identity matrix.
%
%  H = PUPolyMatRand(M,L,N) generates a random paraunitary matrix with a seed
%  value of N. Repeated calls with the same N will results in the same para-
%  unitary matrices.
%
%  H = PUPolyMatRand(M,L,N,'real') generates a real-valyed random paraunitary 
%  while matrix for H = PUPolyMatRand(M,L,N,'complex'), the output will be
%  complex-valued.
%
%  Input parameters:
%     M       spatial dimension of paraunitary matrix
%     L       polynomial order of paraunitary matrix
%     N       seed value for random number generator
%             (optional); default is random;
%     mode    real- or complex valued operation (default is real)
%
%  Output parameter:
%     H       MxMx(L+1) paraunitary matrix
%
%  Reference:
%  [1] P.P. Vaidyanathan: "Multirate Systems and Filter Banks" Prentice
%      Hall, 1993.

%  S. Weiss, University of Strathclyde, 14/12/14
%       amended 31/3/23

%-----------------------------------------------------
%   check for optional seed value
%-----------------------------------------------------
CmplxVld=0;
if nargin>2,
   if nargin>3,
     randn('seed',N);
   else
     if ischar(N)
      mode = N;
     else
      randn('seed',N);
     end
   end;
end;
if strcmp(mode,'complex')==1,
   CmplxVld=1;
end;
      
%-----------------------------------------------------
%   generate L-th order paraunitary matrix
%-----------------------------------------------------
if CmplxVld==1,
   H=randn(M,M) + sqrt(-1)*randn(M,M);
else,   
   H=randn(M,M);
end;   
[u,s,v] = svd(H);
H = u*v';
for i = 1:L,
   % generate elementary paraunitary matrix
   if CmplxVld==1,
      u = randn(M,2)*[1; sqrt(-1)];
   else
      u = randn(M,1);
   end;
   u = u/norm(u,2);
   U(:,:,1)=eye(M)-u*u';
   U(:,:,2)=u*u';
   % apply ith elementary PU matrix
   H = PolyMatConv(H,U);
 end;

