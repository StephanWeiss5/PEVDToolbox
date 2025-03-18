function PHPolyMatInvDemo(sw);
%PHPolyMatInvDemo(SourceSwitch);
%
%  PHPolyMatInvDemo(SourceSwitch) demonstrates the inversion of a 
%  parahermitian matrix using the function PHPolyMatInv(). The input
%  parameter SourceSwitch selects between different sample parahermitian
%  matrices
%
%  Input parameter
%     SourceSwitch    0: a random complex valued 3x3x5 matrix (default) 
%                     1: a simple real valued 3x3x5 matrix
%                     2: a real valued 4x4x5 matrix
%                     3: a simple complex valued 3x3x3 matrix

%  S. Weiss, 6/1/2016

%-----------------------------------------------------------
%  create/pick a parahermitian matrix
%-----------------------------------------------------------
if nargin == 0,         % default for switch selection
   sw=0;
end;
if sw == 0,             % random complex covariance matrix
   A = randn(3,3,3) + sqrt(-1)*randn(3,3,3);;
   R = MIMOConv(A,ParaHerm(A));
elseif sw ==1,          % cov. matrix with 2 off-diag. terms
   R(:,:,1) = zeros(3,3,1);
   R(:,:,2) = zeros(3,3,1);
   R(:,:,3) = diag(ones(3,1));
   R(:,:,4) = zeros(3,3,1);
   R(:,:,5) = zeros(3,3,1);
   R(3,2,1) = .5;
   R(2,3,5) = .5;
   R(1,2,2) = -.4;
   R(2,1,4) = -.4;
elseif sw == 2,
   R(:,:,1) = zeros(3,3,1);
   R(:,:,2) = zeros(3,3,1);
   R(:,:,3) = diag([1 1 3]);
   R(:,:,4) = zeros(3,3,1);
   R(:,:,5) = zeros(3,3,1);
   R(1,2,1) = .5;
   R(2,1,5) = .5;
   R(1,2,2) = -.4;
   R(2,1,4) = -.4;
else                    % cov. matrix with single off-diag. term
   R(:,:,1) = zeros(3,3,1);
   R(:,:,2) = diag(ones(3,1));
   R(:,:,3) = zeros(3,3,1);
   R(3,1,1) = sqrt(-1)*.5;
   R(1,3,3) = -sqrt(-1)*.5;
end;  

S = PHPolyMatInv(R,201);
PolyMatDisplay(abs(PolyMatConv(S,R)));
