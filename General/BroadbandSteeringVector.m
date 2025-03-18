function S = BroadbandSteeringVector(M,alpha,T);
%S = BroadbandSteeringVector(M,alpha,T);
%  
%   S=BroadbandSteeringVector(M,alpha) simulates a broadband steering vector 
%   for a linear equispaced array with critical sampling in both space and time.
%   The number of array elements is M, and the implementation of the broadband
%   steering vector is based on a windowed sinc function [1,2] of order 50.
%    
%   S=BroadbandSteeringVector(M,alpha,T) changes the default value for the order 
%   of the fractional delay filter to 2T. 
%
%   Input parameters:
%      M      number of array elements
%      alpha  angle of incident measured against broadside /[rad]
%      T      half order of fractional delay filter (optional)
%             default  value: T=25;
%
%   Output parameter:
%      S      broadband steering vector (M x 1 x (2T+1) )
%   
%   References:
%   [1]  T.I. Laakso, V. Valimaki, M. Karjalainen, and U.K. Laine: "Splitting 
%        the unit delay," IEEE Signal Processing Magazine, 13(1):30-60, Jan.
%        1996.
%   [2]  J.Selva, "An efficient structure for the design of variable fractional 
%        delay filters based on the windowing method," IEEE Transactions on 
%        Signal Processing, 56(8):3770--3775, Aug. 2008.

%   M. Alrmah, S. Weiss, University of Strathclyde, 24/11/2014
%   updated for delays to be centred w.r.t. the array, S. Weiss, 16/11/22

% input option
if nargin==2,
   T = 25;
end;

% fractional delay between adjacent elements
dT = sin(alpha);

% centre the delay profile
if dT >=0,
   DelayProfile = ((0:(M-1))-(M-1)/2)*dT;
else;
   DelayProfile = -((0:(M-1))-(M-1)/2)*dT;
end;

% maximum shift and effective support of indiviual fractional delay filters
MaxShift = floor(DelayProfile(M)+0.5);
Tr = T-MaxShift;

S = zeros(M,1,2*T+1);
for m =1:M,
   t = (-Tr:Tr);
   FracDelay = DelayProfile(m);
   FracFracDelay = mod(FracDelay+0.5,1)-0.5;
   IntFracDelay = floor(FracDelay+0.5);
   w_Hann = 1 + cos(2*pi*(t-FracFracDelay)/2/(Tr+1));   
   % windowed sinc function
   dummy = sinc(t-FracFracDelay).*w_Hann/sqrt(M);
   % insertion into broadband steering vector
   Shift = MaxShift+IntFracDelay;
   if dT >= 0,
     S(m,1,Shift+(1:2*Tr+1)) = dummy;
   else; 
     S(M-m+1,1,Shift+(1:2*Tr+1)) = dummy;
   end; 
end; 

