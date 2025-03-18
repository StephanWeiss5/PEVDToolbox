function [Pss,Ps,AA] = PssMusic(Q,L,K,angles,epsilon);
% function [Pss,Ps] = PssMusic(Q,L,K,angles);
%
% Polynomial MUSIC algorithm.
%
% Input parameter:
%   Q      polynomial (spatio-temporal) modal matrix
%   L      number of sources
%   K      number of frequency bins evaluated
%   angles possible angles of arrival to be evaluated
%          (measured in rad against broadside)
%   epsilon   regularisation constant to avoid division by zero
%
% Output parameters:
%   Pss    polynomial spatio-spectral music
%   Ps     polynomial spatial music
%
% S. Weiss, University of Strathclyde, 30/1/2013

M = size(Q,1);

Qn = Q(L+1:M,:,:);
Rn = PolyMatConv(ParaHerm(Qn),Qn);


for i = 1:size(angles,2),
    X = BroadbandSteeringVector(M,angles(i),100);
    Y =  PolyMatConv(ParaHerm(X),PolyMatConv(Rn,X));
    dummy2 = shiftdim(Y,1);
    AA(i,:) = dummy2;
end;

if K < size(AA,2),
   disp(['warning: FFT length is shorter than auto-correlation sequence ' ...
         'in function PssMusic()']);
end;

Pss = 1./(abs(fft(AA,K,2).') + epsilon);
Ps = 1./abs(AA(:,(size(AA,2)+1)/2) + epsilon);

