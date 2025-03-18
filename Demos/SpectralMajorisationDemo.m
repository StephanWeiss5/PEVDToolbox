function SpectralMajorisationDemo();
%SpectralMajorisationDemo()
%
%  This demo demonstrated spectral majorisation achieved by a number of 
%  iterative PEVD algorithms after 50 and 200 iterations for a situation 
%  with known ground truth.
%  The ground truth is provided by a source model generated by the 
%  function GenerateMIMOSources(), which here contains 8 sources with 
%  spectrally majorised PSDs and a paraunitary convolutive mixing matrix
%  of order 16.
%
%  The results reflect approximately those shown in Figures 4 and 5 of [1].
%
%  Reference:
%  
%  [1] S. Redif, S. Weiss, J.G. McWhirter: "Sequential Matrix Diagonalisation 
%      Algorithms for Polynomial EVD of Parahermitian Matrices", to appear in
%      IEEE Transactions on Signal Processing.

% S. Weiss, University of Strathclyde, 2/11/2014

%-----------------------------------------------------
%  Ground truth and simulation parameters
%-----------------------------------------------------
% set seeds for a specific source model
randn('seed',10); rand('seed',10);
% generate source model with the following parameters:
L = 8;         %  # of sources
P = 16;        %  order of source innovation filter
M = 8;         %  # of sensors
K = 16;        %  order of paraunitary mixing matrix
gamma = 0.08;  %  max radii of zeros
% convolutive mixing model
[H,D,F] = GenerateMIMOSources(L,P,M,K,gamma);
% PSDs of ground truth sources
P = PolyMatDiagSpec(D,1024);    
% space-time covariance matrix at sensors
R = PolyMatConv(H,PolyMatConv(D,ParaHerm(H)));
% some display parameters
LWidth = 1;    % line widths for spectra
LWidth2 = 3;   % line width for ground truth
ShadeColour=[.7 .7 .7];  % grey scale for ground truth
w= (0:1023)/512;   % frequency scale

%-----------------------------------------------------
%   SMD, ME-SMD, SBR2 and SBR2C after 50 iterations
%-----------------------------------------------------
figure(1);
clf;
[H_est,D_est] = SMD(R,50,0.00000001,0.00000001,'SMD');
P_est = PolyMatDiagSpec(D_est,1024);
subplot(221);
% entries for legend
plot([-2 -1],[1 2],'b-o','LineWidth',LWidth);
hold on;
plot([-2 -1],[1 2],'r-.*','LineWidth',LWidth);
plot([-2 -1],[1 2],'k--s','LineWidth',LWidth);
plot([-2 -1],[1 2],'-+','Color',[0 .5 0],'LineWidth',LWidth);
AlgorithmPSDs50Plot(P_est,P,LWidth,LWidth2,ShadeColour);
% legend and labels
dummy = legend('$l=1$','$l=2$','$l=3$','$l=4$');
set(dummy,'interpreter','latex','location','SouthWest');
text(0.05,1,'(a)');
title('SMD after 50 iterations');

% ME-SMD
subplot(222);
[H_est,D_est] = SMD(R,50,0.00000001,0.00000001,'MESMD');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs50Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(b)');
title('ME-SMD after 50 iterations');

% SBR2
subplot(223);
[H_est,D_est] = SBR2(R,50,0.00000001,0.00000001,'SBR2');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs50Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(c)');
title('SBR2 after 50 iterations');

% SBR2C
subplot(224);
[H_est,D_est] = SBR2(R,50,0.00000001,0.00000001,'SBR2C');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs50Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(d)');
title('SBR2C after 50 iterations');

%-----------------------------------------------------
%   SMD, ME-SMD, SBR2 and SBR2C after 200 iterations
%-----------------------------------------------------
figure(2);
clf;
[H_est,D_est] = SMD(R,200,0.00000001,0.00000001,'SMD');
P_est = PolyMatDiagSpec(D_est,1024);
subplot(221);
% entries for legend
plot([-2 -1],[1 2],'b-o','LineWidth',LWidth);
hold on;
plot([-2 -1],[1 2],'r-.*','LineWidth',LWidth);
plot([-2 -1],[1 2],'k--s','LineWidth',LWidth);
plot([-2 -1],[1 2],'-+','Color',[0 .5 0],'LineWidth',LWidth);
plot([-2 -1],[1 2],'b-.d','LineWidth',LWidth);
plot([-2 -1],[1 2],'r--x','LineWidth',LWidth);
plot([-2 -1],[1 2],'k-o','LineWidth',LWidth);
plot([-2 -1],[1 2],'-.*','Color',[0 .5 0],'LineWidth',LWidth);
AlgorithmPSDs200Plot(P_est,P,LWidth,LWidth2,ShadeColour);
dummy = legend('$l=1$','$l=2$','$l=3$','$l=4$','$l=5$','$l=6$','$l=7$','$l=8$');
set(dummy,'interpreter','latex','location','SouthWest');
text(0.05,1,'(a)');
title('SMD after 200 iterations');

% ME-SMD
subplot(222);
[H_est,D_est] = SMD(R,200,0.00000001,0.00000001,'MESMD');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs200Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(b)');
title('ME-SMD after 200 iterations');

% SBR2
subplot(223);
[H_est,D_est] = SBR2(R,200,0.00000001,0.00000001,'SBR2');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs200Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(c)');
title('SBR2 after 200 iterations');

% SBR2C
subplot(224);
[H_est,D_est] = SBR2(R,200,0.00000001,0.00000001,'SBR2C');
P_est = PolyMatDiagSpec(D_est,1024);
AlgorithmPSDs200Plot(P_est,P,LWidth,LWidth2,ShadeColour);
text(0.05,1,'(d)');
title('SBR2C after 200 iterations');


%=============================================
function AlgorithmPSDs50Plot(P_est,P,LWidth,LWidth2,ShadeColour);
% Plot power spectral densities of first four columns in P and P_est
% background
w = (0:1023)/512;
plot(w,10*log10(abs(P(:,1))),'Color',ShadeColour,'LineStyle','-','LineWidth',LWidth2);
hold on;
plot(w,10*log10(abs(P(:,2))),'Color',ShadeColour,'LineStyle','-.','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,3))),'Color',ShadeColour,'LineStyle','--','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,4))),'Color',ShadeColour,'LineStyle','-','LineWidth',LWidth2);
% lines
plot(w,10*log10(abs(P_est(:,1))),'b-','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,2))),'r-.','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,3))),'k--','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,4))),'-','Color',[0 .5 0],'LineWidth',LWidth);
% markers
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,1))),'bo','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,2))),'r*','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,3))),'ks','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,4))),'+','Color',[0 .5 0],'LineWidth',LWidth);
% legend and labels
xlabel('norm. angular frequency $\Omega/\pi$','fontsize',14,...
  'Interpreter','Latex');
ylabel('$10\log_{10}S^{(50)}_{s,l}(e^{j\Omega})$ / [dB]','fontsize',...
  14,'Interpreter','Latex');
axis([0 2 -8 2]);
grid on;

%=============================================
function AlgorithmPSDs200Plot(P_est,P,LWidth,LWidth2,ShadeColour);
% Plot power spectral densities of all eight columns in P and P_est
% background
w= (0:1023)/512;

% background
plot(w,10*log10(abs(P(:,1))),'Color',ShadeColour,'LineStyle','-','LineWidth',LWidth2);
hold on;
plot(w,10*log10(abs(P(:,2))),'Color',ShadeColour,'LineStyle','-.','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,3))),'Color',ShadeColour,'LineStyle','--','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,4))),'Color',ShadeColour,'LineStyle','-','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,5))),'Color',ShadeColour,'LineStyle','-.','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,6))),'Color',ShadeColour,'LineStyle','--','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,7))),'Color',ShadeColour,'LineStyle','-','LineWidth',LWidth2);
plot(w,10*log10(abs(P(:,8))),'Color',ShadeColour,'LineStyle','-.','LineWidth',LWidth2);
% lines
plot(w,10*log10(abs(P_est(:,1))),'b-','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,2))),'r-.','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,3))),'k--','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,4))),'-','Color',[0 .5 0],'LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,5))),'b-.','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,6))),'r--','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,7))),'k-','LineWidth',LWidth);
plot(w,10*log10(abs(P_est(:,8))),'-.','Color',[0 .5 0],'LineWidth',LWidth);
% markers
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,1))),'bo','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,2))),'r*','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,3))),'ks','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,4))),'+','Color',[0 .5 0],'LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,5))),'bd','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,6))),'rx','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,7))),'ko','LineWidth',LWidth);
plot(w(64:128:end),10*log10(abs(P_est(64:128:end,8))),'g*','Color',[0 .5 0],'LineWidth',LWidth);
% axes
xlabel('norm. angular frequency $\Omega/\pi$','fontsize',14,...
  'Interpreter','Latex');
ylabel('$10\log_{10}S^{(200)}_{s,l}(e^{j\Omega})$ / [dB]','fontsize',...
  14,'Interpreter','Latex');
axis([0 2 -14 2]);
grid on;

