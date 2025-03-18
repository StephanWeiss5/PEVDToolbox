function PolyMatDisplay(H,t,mode);
%PolyMatDisplay(H,t,mode);
%  
%   Display the MIMO system or polynomial matrix represented by H. If H is
%   of dimension K x M x L, then Matlab will produce KxM subplots, each con-
%   taining an FIR response of length L.
%
%   In general, the system H(z) represented by H will be of the form
%      H(z) = H_{-L1} z^{L1} + ... + H_{-1} z + H0 + 
%                + H1 z^{-1} + ... + H_{L2} z^{-L2}
%   such that the total support length is L = L1+L2+1.
%
%   PolyMatDisplay(H) will plot the FIR filters along an x-axis range [0;L-1].
%   PolyMatDisplay(H,t) will plot the FIR filters against an axis t, whereby the
%   length of t must equal L.
% 
%   The input H is assumed to be real valued. If complex valued, only the 
%   absolute values of the coefficients will be plotted.
%
%   The function produces a Matlab figure with KxM subplots each containing
%   a stem plot. The axes are returned unlabelled.
%
%   PolyMatDisplay(H,t,'hold') plots H over t in red * markers over an existing
%   plot.
%
%   Input parameters
%      H      KxMXL MIMO system matrix or polynomial matrix
%      t      L-dimensional axis vector (optional)
%             default:  t=(0:(L-1));
%      mode   'hold' to plot over an existing plot (optiona)
%              default if off
%
%   Output parameters: none

%    S Weiss, Univ. of Strathclyde, 30/7/14
%        added 'hold' option 21/12/23
  
hold off;  
[K,L,M] = size(H);

%------------------------------------------------------------------------------
% check dimensions and options
%------------------------------------------------------------------------------
if nargin==1,
  t = (0:M-1);
end;  
if nargin<=2,
   hold off;     
   HoldMode = 0;
else
   if mode=='hold',
       HoldMode = 1;
   else
       error('plot option undefined');    
   end;    
end;        
if length(t)~=M,
   error('dimension mismatch for input parameters to function MIMODisplay()');
end;

%------------------------------------------------------------------------------
% check if input matrix is real valued
%------------------------------------------------------------------------------
if isreal(H)~=1,
   H = abs(H);    % if not, plot absolute value of coefficients
end;

Hmin = min(min(min(H)));
Hmax = max(max(max(H)));
if Hmin>0, 
  Hmin=0; 
end;
i = 0;
for k = 1:K,              % rows
  for l = 1:L,            % columns
    i = i + 1;
    subplot(K,L,i);
    if HoldMode==0,
       stem(t,shiftdim(H(k,l,:),1));
       axis([t(1)-.1 t(end)+.1 Hmin*1.1 Hmax*1.1]);
    else
       hold on; 
       Aold = axis;
       Anew =   [t(1)-.1 t(end)+.1 Hmin*1.1 Hmax*1.1];
       A = [min(Aold(1),Anew(1)) max(Aold(2),Anew(2)) min(Aold(3),Anew(3)) max(Aold(4),Anew(4))];
       plot(t,shiftdim(H(k,l,:),1),'r*');
       axis(A);
    end;   
  end;
end;
