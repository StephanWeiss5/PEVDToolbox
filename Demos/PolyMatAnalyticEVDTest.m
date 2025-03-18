% PolyMatAnalyticEVDTest.m
%
% Demonstrate the use of the functions PolyMatAnalyticEigValues() and 
% PolyMatAnalyticEigVectors().
%
% S. Weiss 11/5/2022

clear all; close all;

%--------------------------------------------------------------
%   parameters
%--------------------------------------------------------------
M = 3;       % spatial dimension
Ll = 3;      % length of eigenvalues (must be odd)
Lq = 3;      % length of eigenvectors

%--------------------------------------------------------------
%   construct a parahermitian matrix
%--------------------------------------------------------------
Lambda=zeros(M,M,Ll); 
for i = 1:M,
  dummy = randn((Ll+1)/2,1);
  Lambda(i,i,:) = conv(dummy,conj(flipud(dummy)));
end;   
Q = zeros(3,3,1);
Q(:,:,1) = eye(3);
for i = 1:Lq-1,
   dummy = randn(M,1);
   dummy = dummy/norm(dummy);
   U = zeros(M,M,2);
   U(:,:,1) = eye(3)-dummy*dummy';
   U(:,:,2) = dummy*dummy';
   Q = PolyMatConv(Q,U);
end;   
R = PolyMatConv(Q,PolyMatConv(Lambda,ParaHerm(Q)));

%--------------------------------------------------------------
%   Analytic EVD
%--------------------------------------------------------------
L_hat = PolyMatAnalyticEigValues(R);
Q_hat = PolyMatAnalyticEigVectors(R,L_hat);
Lambda_hat = zeros(M,M,size(L_hat,2));
for m = 1:M,
   Lambda_hat(m,m,:) = L_hat(m,:);
end;
   
%--------------------------------------------------------------
%   Diagnostics
%--------------------------------------------------------------
figure(1); PolyMatDisplay(real(Lambda));
figure(2); PolyMatDisplay(real(Q));
disp('ground truth in Figs. 1 and 2');
figure(3); PolyMatDisplay(real(Lambda_hat));
figure(4); PolyMatDisplay(real(Q_hat));
disp('extractions in Figs. 3 and 4');

