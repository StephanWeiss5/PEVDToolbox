function ConvergenceDemo(SourceSwitch);
%ConvergenceDemo(SourceSwitch);
%
%  ConvergenceDemo(SourceSwitch) compares the polynomial EVD approximations by 
%  the functions SBR2() and SMD() for different types of parahermitian matrices. 
%  The demo provides the reduction in off-diagonal energy in increments of 5 
%  iteration steps.
% 
%  Theinput parameter SourceSwitch can be selected as:
%   
%  ComparisonDemo('Random') uses an arbitrary complex valued 5x5 polynomial 
%  parahermitian matrix of order 10.
%
%  ComparisonDemo('Simple1') compares decompositions of a 3x3 real valued 
%  polynomial parahermitian matrix R(z) of order 4:
%                | 1     -.4*z       0     |
%       R(z) =   | -.4*z   1     .5*z^{-2} |
%                | 0     .5*z^2      3     |     
%
%  ComparisonDemo('Simple2') compares the decompositions of a 5x5 complex valued 
%  polynomial parahermitian matrix R(z) of order 4:
%                |    1        .5*z^2  -.4j*z^{-1}    0   .2j*z^2 |
%                | .5*z^{-2}      1         0         0       0   |
%       R(z) =   |   .4j*z        0         3      .3z^{-1}   0   |
%                |    0           0       .3*z       .5       0   |
%                |-.2j*z^{-2}     0         0         0      .25  |
%  This matrix admits a PEVD with a non-polynomial Gamma(z).
%
%  The function generates a plot showing the remaining off-diagonal power, 
%  normalised w.r.t. the total power of the matrix; this indicates how well the
%  two algorithms perform their diagonalisation task. A second figure compares 
%  the off-diagonal power in relation to the order of the paraunitary matrices 
%  required in order to accomplish this decomposition. 
%
%  Input parameter:
%       SourceSwitch       selects parahermitian matrx to be decomposed
%                          default: 'Random'
%  Output parameter: none

% S. Weiss, University of Southampton, 8/10/2014

if nargin==0,
   SourceSwitch='Random';
end;

disp('SBR2 and SMD Convergence Comparison Demo');
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
%  Set incremental parameters 
%---------------------------------------
% parameters
maxiter = 5;                  % max number of iterations per update
blocks  = 20;
epsilon = 10^(-10);           % alternative stopping criterion
                              %    ensure that this is not invoked
mu = 0;                       % truncate true zeroes
GammaSMD = R;                 % initialisations
GammaSBR2 = R;
N1 = PolyMatNorm(R);
N2 = PolyMatNorm(R,'OffDiag');
OffDiagNormSMD(1) = N2/N1;    % normalised remaining off-diag. energy
OffDiagNormSBR2(1) = N2/N1;        
OrderSMD(1)=0;
OrderSBR2(1)=0;

%---------------------------------------
%  Recursively decompose in blocks of 5 iterations
%---------------------------------------
for n = 2:blocks,
   % SMD iteration
   [H,GammaSMD] = SMD(GammaSMD,maxiter,epsilon,mu,'SMD');
   OffDiagNormSMD(n) = PolyMatNorm(GammaSMD,'OffDiag')/N1;
   OrderSMD(n) = OrderSMD(n-1)+size(H,3)-1;
   % SBR2 iteration
   [H,GammaSBR2] = SBR2(GammaSBR2,maxiter,epsilon,mu,'SBR2');
   OffDiagNormSBR2(n) = PolyMatNorm(GammaSBR2,'OffDiag')/N1;
   OrderSBR2(n) = OrderSBR2(n-1)+size(H,3)-1;
end;

%---------------------------------------
%  Display resulting matrices 
%---------------------------------------
% display normalised remaining off-diagonal power vs iterations 
figure(1); clf; 
plot(1:maxiter:maxiter*blocks,5*log10(OffDiagNormSBR2),'bo');
hold on;
plot(1:maxiter:maxiter*blocks,5*log10(OffDiagNormSMD),'r*');
legend('SBR2','SMD');
xlabel('iterations');
ylabel('norm. remaining off-diagonal power / [dB]');
disp('Figure 1: norm. remaining off-diagonal power vs iterations');

% display normalised remaining off-diagonal power vs paraunitary order
figure(2); clf;
plot(OrderSBR2,5*log10(OffDiagNormSBR2),'bo');
hold on;
plot(OrderSMD,5*log10(OffDiagNormSMD),'r*');
legend('SBR2','SMD');
xlabel('paraunitary order (without truncation)');
ylabel('norm. remaining off-diagonal power / [dB]');
disp('Figure 2: norm. remaining off-diagonal power vs paraunitary order');
