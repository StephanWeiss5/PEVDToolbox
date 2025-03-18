function [q,chi] = PhaseSmoothing(q_in,P,q_opt);
% [q,chi] = PhaseSmoothing(q_in,P);
%
%   With q_in an MxNfft matrix, this function aims to change the phase of every
%   column of q_in such that its inverse DFT yields a support that is as 
%   compact as possible. This is equivalent to smoothing the phase across the 
%   column vectors of q_in. 
%
%   This smoothing is performed by minimising the power in the P-the derivative 
%   of the Dirichlet interpolations through all M functions.   
%
%   Input parameters:
%      q_in       MxNfft matrix, representing a frequency-dependent M-element 
%                    vector
%      P          derivative order for power calculation
%      q_opt      optimal phase adjustment (optional & for comparison only)
%
%   Output parameters:
%      q          frequency-dependent vector with phase adjustment, such that 
%                 ifft(q,[],2) has maximally short support.
%      chi        smoothness metric at the solution

%   S. Weiss, UoS, 22/7/2020
%         updated with modulation correction, 31/7/2020
%         updated to account for correct Hessian, and working with only (Nfft-1) angles, 5/1/22

%-----------------------------------------
%  derivative and internal parameters
%-----------------------------------------
DISPLAY = 0;        % plot intermediate results
if DISPLAY==1, disp('***** in routine PhaseSmoothing() *****'); end;
[M,Nfft] = size(q_in);

%-----------------------------------------
%  set up smoothness metric
%-----------------------------------------
c = (0:(Nfft-1)).^(P);
W = dftmtx(Nfft)/Nfft;
C = W*diag(c.^2)*W'; 

%-----------------------------------------
%  display interpolations of starting point
%-----------------------------------------
if DISPLAY==1,
   figure(1); clf;
   for m = 1:M,
     f1 = fft(ifft(q_in(m,:)),10*Nfft);
     f =  (0:(Nfft-1))/Nfft;
     f10 = (0:(10*Nfft-1))/10/Nfft;
     plot3(f10,real(f1),imag(f1)); 
     if m ==1, hold on; end;
     plot3(f,real(q_in(m,:)),imag(q_in(m,:)),'*');
   end;

   if nargin>2,
     figure(2); clf;
     for m = 1:M,
       f1 = fft(ifft(q_opt(m,:)),10*Nfft);
       f =  (0:(Nfft-1))/Nfft
       f10 = (0:(10*Nfft-1))/10/Nfft;
       plot3(f10,real(f1),imag(f1)); 
       if m ==1, hold on; end;
       plot3(f,real(q_opt(m,:)),imag(q_opt(m,:)),'*');
     end;
   end;
end;

%-----------------------------------------
%  starting metrics
%-----------------------------------------
if nargin > 2,
  chi_opt = norm(ifft(q_opt,Nfft,2)*diag(c),'fro').^2;
else
  chi_opt = 1.0;
end;

%-----------------------------------------
%  initialisation
%-----------------------------------------
q_in = q_in*exp(-j*angle(q_in(1,1)));   % first reference element is real
q = zeros(M,Nfft);

%-----------------------------------------
%  iteration
%-----------------------------------------
q = q_in;
prof = zeros(Nfft,1);
chi(1) = norm(ifft(q,Nfft,2)*diag(c),'fro').^2/chi_opt;
if DISPLAY == 1, disp(sprintf('iteration 0, relative smoothness  %f',chi(1))); end;
D = zeros(Nfft,Nfft);
for m = 1:M, D = D + diag(q(m,:))'*C*diag(q(m,:)); end;
phi = zeros(Nfft,1);
mu = 2; crit1 = 1; crit2 = 1; n = 1; chi_old = 1/eps; Mod = 0;
while crit1 == 1, 
  n = n + 1;
  % update as proposed in SSPD'20 paper
  Phi = diag(exp(j*phi));
  a = exp(j*phi);
  Ar = [zeros(1,Nfft-1); diag(a(2:Nfft))];
  dummy = Ar'*D*a;
  H = 2*real(Ar'*D*Ar); %- 2*diag(real(dummy));
  Hinv = pinv(H);
%  phi = phi - mu*Hinv*imag(conj(Phi)*D*Phi*ones(Nfft,1));
  phi(2:Nfft) = phi(2:Nfft) - mu*Hinv*imag(dummy);
%    phi = phi - mu*imag(conj(Phi)*D*Phi*ones(Nfft,1));
  phi = phi - phi(1);    % first reference element is real
  for k = 1:Nfft,
    q(:,k) = q_in(:,k)*exp(j*phi(k));
  end;
  % the two following evaluations of the cost are equivalent
  %  chi(n) = ones(1,Nfft)*conj(Phi)*D*Phi*ones(Nfft,1);   
  chi(n) = norm(ifft(q,Nfft,2)*diag(c),'fro').^2/chi_opt;
  ratio = chi(n)/chi(n-1);
  if DISPLAY==1,
    disp(sprintf('iteration %d, relative smoothness  %f',[(n-1) chi(n)]));
    disp(sprintf('ratio %f',ratio));
    f2 = (0:Nfft*10-1)/10/Nfft;
    Q1 = fft(ifft(q,Nfft,2),10*Nfft,2);
    figure(1); clf;
    for m = 1:M,
      plot3(f2,real(Q1(m,:)),imag(Q1(m,:)),'b-');
      if m == 1, hold on; end; 
      plot3((0:Nfft-1)/Nfft,real(q(m,:)),imag(q(m,:)),'o'); 
    end;
    drawnow;
    prof = sum(abs(ifft(q,Nfft,2)),1);
    figure(3); stem(1:Nfft,prof);
    drawnow; pause(0.25);
  end;
  % now need to check whether we have ended up in a critical point, and require modulation
  if (crit2 == 1) & (n<Nfft+10),
    if (ratio<1) & (ratio>0.99)   % learning curve has levelled out
      if (chi(n) < chi_old),   % smaller than at previous critical point -> modulation
        phi_old = phi;
        chi_old = chi(n);
        phi = phi+(2*pi/Nfft*(0:(Nfft-1))');
        if DISPLAY==1, disp('modulated'); end;
        Mod = 1;
      end;
    elseif (ratio>1000)&(Mod==1),
      % cost function is now worse than at the previous critical point
      %   [NEED TO FIND A MORE RELIABLE TEST FOR THE RATIO!]
        phi = phi_old;  % undo last modulation
        crit2 = 0;      % suppress further modulation steps
        if DISPLAY == 1, disp('undo last modulation step'); end;
    else
      % do nothing
    end;
  else  % now iterate a bit longer until desired accuracy is reached
    if DISPLAY == 1, disp('final improvements'); end;
    if (ratio>0.99999), crit1 = 0; end;     % max number of iterations
  end;
end;
for k = 1:Nfft,
  q(:,k) = q_in(:,k)*exp(j*phi(k));
end;
