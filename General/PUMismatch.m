function e = PUMismatch(Q);
% PUMismatch(Q)
%
% e = PUMismatch(Q) measures the mismatch of the polynomial Q to a paraunitary
% system by evaluating Q(z) Q^P(z), and measuring the least squares mismatch 
% to an identity matrix. The least squares error is returned in the variable e.
%
% Input parameter
%    Q         square polynomial matrix (MxMxL)
%
% Output parameter
%    e         paraunitarity error

% S. Weiss, UoS, 20/1/2025

 [M,~,LQ] = size(Q);
 Nfft = 2^(ceil(log2(LQ))+1);
 Qf = fft(Q,Nfft,3);
 e = 0;
 for k = 1:Nfft,
    e = e + norm(eye(M) - Qf(:,:,k)*Qf(:,:,k)','fro')^2; 
 end;   
 e = e/Nfft;
