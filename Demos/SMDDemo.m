function SMDDemo(SourceSwitch);
%SMDDemo(SourceSwitch);
%
%  SMDDemo(SourceSwitch) demonstrates the polynomial EVD approximated by 
%  the function SMD() for different type of parahermitian matrices. The
%  input parameter SourceSwitch can be selected as:
%
%  SMDDemo('Random') produces an abritrary complex valued 5x5 polynomial 
%  parahermitian matrix of order 10.
%
%  SMDDemo('Simple1') decomposes a 3x3 real valued polynomial parahermitian 
%  matrix R(z) of order 4:
%                | 1     -.4*z       0     |
%       R(z) =   | -.4*z   1     .5*z^{-2} |
%                | 0     .5*z^2      3     |     
%
%  SMDDemo('Simple2') decomposes a 5x5 complex valued polynomial 
%  parahermitian matrix R(z) of order 4:
%                |    1        .5*z^2  -.4j*z^{-1}    0   .2j*z^2 |
%                | .5*z^{-2}      1         0         0       0   |
%       R(z) =   |   .4j*z        0         3      .3z^{-1}   0   |
%                |    0           0       .3*z       .5       0   |
%                |-.2j*z^{-2}     0         0         0      .25  |
%  This matrix admits a PEVD with a non-polynomial Gamma(z).
%
%  The function generates a number of outputs to highlight diagonalisation
%  and spectral majorisation.
%
%  Input parameter:
%       SourceSwitch       selects parahermitian matrx to be decomposed
%                          default: 'Random'
%  Output parameter: none

% S. Weiss, University of Southampton, 3/10/20

if nargin==0,
   SourceSwitch='Random';
end;

disp('Sequential Matrix Diagonalisation (SMD) Demo');
disp('--------------------------------------------');

%---------------------------------------
%  Define Scenario
%---------------------------------------
% create parahermitian matrix to be decomposed
if strcmp(SourceSwitch,'Random')==1,
   A = randn(5,5,6) + sqrt(-1)*randn(5,5,6);;
   R = PolyMatConv(A,ParaHerm(A));
elseif strcmp(SourceSwitch,'Simple1')==1,
   R        = zeros(3,3,5);
   R(:,:,3) = diag([1 1 3]);
   R(1,2,1) = .5;
   R(2,1,5) = .5;
   R(1,2,2) = -.4;
   R(2,1,4) = -.4;
elseif strcmp(SourceSwitch,'Simple2')==1,
   R        = zeros(5,5,5);
   R(:,:,3) = diag([1 1 3 .5 .25]);
   R(1,2,1) = .5;
   R(2,1,5) = .5;
   R(1,3,4) = -.4*sqrt(-1);
   R(3,1,2) = .4*sqrt(-1);
   R(1,5,1) = .2*sqrt(-1);
   R(5,1,5) = -.2*sqrt(-1);
   R(4,3,2) = .3;
   R(3,4,4) = .3;
else                    
   error('option for input parameter SourceSwitch not implemented');
end;  

%---------------------------------------
%  Calculate PEVD via SBR2 
%---------------------------------------
% parameters
maxiter = 100;                % max number of iterations
epsilon = 0.001;              % alternative stopping criterion
mu = 0;                       % truncate true zeroes
% decomposition
[H,Gamma] = SMD(R,maxiter,epsilon,mu,'SMD');

%---------------------------------------
%  Display resulting matrices 
%---------------------------------------
% display original matrix
figure(1); clf; 
L=size(R,3);
t=(-(L-1)/2:(L-1)/2);
PolyMatDisplay(abs(R),t);
disp('Figure 1: original parahermitian matrix R(z)');

% display paraunitary matrix
figure(2); clf;
L=size(H,3);
PolyMatDisplay(abs(H));
disp('Figure 2: paraunitary matrix H(z)');

% display diagonalised matrix
figure(3); clf;
L=size(Gamma,3);
t=(-(L-1)/2:(L-1)/2);
PolyMatDisplay(abs(Gamma),t);
disp('Figure 3: diagonalised matrix Gamma(z)');

% demonstrate approximate spectral majorisation
figure(4); clf;
Ndft = max([256,size(Gamma,3)]); 
P = PolyMatDiagSpec(Gamma,Ndft);
plot((0:Ndft-1)/Ndft,10*log10(abs(P)));
xlabel('normalised angular frequency \Omega/(2\pi)');
ylabel('power spectral densities / [dB]');
disp('Figure 4: PSDs of diagonal elements of Gamma(z)');
   
%---------------------------------------
%  Check result quantitatively 
%---------------------------------------
% remaining off-diagonal energy
N1 = PolyMatNorm(R);
N2 = PolyMatNorm(Gamma,'OffDiag');
disp(sprintf('ratio off-diagonal/total energy: %f (%f dB)',...
     [N2/N1, 10*log10(N2/N1)]));

% accuracy of decomposition and reconstruction
R2 = PolyMatConv(ParaHerm(H),PolyMatConv(Gamma,H));
L  = size(R,3);
L2 = size(R2,3);
Indices = (L2+1)/2 + ( -(L-1)/2:(L-1)/2 );
N3 = PolyMatNorm(R-R2(:,:,Indices));
disp(sprintf('error  |R(z) -  H~(z)Gamma(z)H(z)| = %f (%f dB)',...
     [N3, 10*log10(N3)]));
