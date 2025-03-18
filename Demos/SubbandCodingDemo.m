function SubbandCodingDemo(SourceSwitch)
%SubbandCodingDemo(SourceSwitch)
%
%  Demonstrates the use of SBR2 and SBR2C to the problem of subband coding. 
%  The optional subband coder is given by the paraunitary matrix extracted 
%  by a PEVD. 
%  
%  SubbandCodingDemo(SourceSwitch) provides two demonstrations selected 
%  through the variable SourceSwitch.
%
%  SubbandCodingDemo('MA14') generates a random moving average innovation 
%  filter of order 14, with coefficients drawn from a complex Gaussian 
%  distribution, from which an auto-correlation sequence is derived. This 
%  process is assumed to be multiplexed into 4 subchannels. 
%
%  SubbandCodingDemo('AR4') approximates an autoregressive (RA) process
%  with two complex pole pairs. This AR4 process is approximated by 
%  finite length auto-correlation of order 200.
%
%  SubbandCodingDemo('Simple') uses a short auto-correlation sequence of 
%  order 6, with a demultiplexing into 3 channels.
%
%  The examples 'MA14' and 'AR4' are similar to the examples in [1].
%
%  Reference:
%
%  [1] S. Redif, J.G. McWhirter, and S. Weiss, "Design of FIR Paraunitary 
%      Filter Banks for Subband Coding Using a Polynomial Eignvalue Decompo- 
%      sition," IEEE Transactions on Signal Processing, vol. 59, no. 11, 
%      pp. 5253-5264, Nov 2011.

%  S Weiss, Univ. of Strathclyde, 27/8/14 

if nargin==0,
   SourceSwitch='Simple';
end;

%---------------------------------------
%  Define Scenario
%---------------------------------------
% create autocorrelation sequence for scenatio
if strcmp(SourceSwitch,'MA14')==1,
   a = (randn(14,1)+sqrt(-1)*randn(14,1))/sqrt(2);
   r = conv(a,flipud(conj(a)));
   M = 4;
elseif strcmp(SourceSwitch,'AR4')==1, 
   poles=[0.9*exp(j*0.6283), 0.9*exp(-j*0.6283), 0.85*exp(j*2.8274), 0.85*exp(-j*2.8274)];
   x = [1; zeros(99,1)];
   a = filter(1,poly(poles),x);
   r = conv(a,flipud(conj(a)));
   M = 4;
elseif strcmp(SourceSwitch,'Simple')==1,
   r = [-j/4 .5 -j 2 j .5 j/4];
   M = 3;
else                    
   error('option for input parameter SourceSwitch not implemented');
end;  

% create parahermitian matrix from autocorrelation sequence
R = MuxPolyCovMat(r,M);

%---------------------------------------
%  Perform decompositions
%---------------------------------------
% SBR2
[H,Gamma]=SBR2(R,50,0.00000001,0,'SBR2');
P1 = PolyMatDiagSpec(Gamma,1024);
CG1 = CodingGain(Gamma);
% SBR2C (SBR2 version optimised for coding gain)
[H,Gamma_C]=SBR2(R,50,0.00000001,0,'SBR2C');
CG2 = CodingGain(Gamma_C);
P2 = PolyMatDiagSpec(Gamma_C,1024);

%---------------------------------------
%  Display results
%---------------------------------------
% numerical results
disp('Subband Coding Demo');
disp('-------------------');
disp(sprintf('  coding gain with SBR2:  %f',CG1));
disp(sprintf('  coding gain with SBR2C: %f',CG2));

% graphical output
W=(0:1023)/512;
figure(1); clf;
plot(W,10*log10(abs(P1)));
hold on;
plot(W,10*log10(abs(P2)),'--');
xlabel('normalised angular frequency \Omega/\pi');
ylabel('power spectral densities / [dB]');
disp('power spectral densities are displayed in Figure 1:')
disp('  solid lines: PSDs achieved by SBR2'); 
disp('  dashed lines: PSDs achieved by SBR2C'); 

