function [xi_1,xi_2,VarE] = SpaceTimeCovMatMSE(R,N,T,VecMode);
%SpaceTimeCovMatMSE(R,N,T,mode);
%
%  SpaceTimeCovMatMSE(R,N,T) determines the mean square estimation error of 
%  a space time covariance matrix R over a support length of (-T ... +T) 
%  given a data set containing N snap shots, based on the analysis in [1].
%
%  Via [xi_1,xi_2] = SpaceTimeCovMatMSE(R,N,T), this function returns two
%  arguments: the estimation error due to the finite sample set in xi_1, 
%  and a truncation error in xi_2. If 2T+1 exceeds the support of R, then 
%  xi_2 will be zero.
%
%  [xi_1,xi_2,VarE] = SpaceTimeCovMatMSE(R,N,T) additionally returns the
%  variance of the estimation error for each element of the estimate.
%     
%  [xi_1,xi_2,VarE] = SpaceTimeCovMatMSE(R,N,T,'vec') returns vectors of 
%  dimension T for xi_1 and xi_2, containing errors for lags tau, 0<tau<=T.
%
%  Input parameters:
%    R         MxMxL space-time covariance matrix
%    N         number of snapshots
%    T         maximum lag range over which estimate is calculated
%    mode      optional, 'vec' for vectorial rather than total error terms  
%     
%  Output parameters
%    xi_1      estimation error due to finite sample size N
%    xi_2      truncation error, xi_2>0 in case T<L
%    VarE      MxMx(2T+1) variance of estimated space time covariance
%                matrix, element by element
%
%  Reference:
%    [1] C. Delaosa, J. Pestana, N. Goddard, S. Somasundaram, S. Weiss:
%        "Sample Space-Time Covariance Matrix Estimation", submitted to 
%        IEEE ICASSP, Brighton, UK, May 2019.
  
% S. Weiss, UoS, 28/10/18
%    return a vector of errors       06/05/2022

%------------------------------------------------------------------------------
%  parameters
%------------------------------------------------------------------------------
[M,~,L] = size(R);
Tmax = (L-1)/2;            % lag range of R is (-Tmax ... +Tmax)
ComplexMode = isreal(R);   % 0: complex;  1: real
VectorMode=0;
if nargin>3,
   if VecMode=='vec',
      VectorMode=1;
   end;
end;      
   
%------------------------------------------------------------------------------
%  check for truncation
%------------------------------------------------------------------------------
if VectorMode==0,
  xi_2 = 0.0;
  if T<Tmax,                 % only evaluate positive lag range
    %   disp('truncation error incurred');
    for t = T+1:Tmax,   
      xi_2 = xi_2 + 2*norm(R(:,:,Tmax+1+t),'fro')^2;
    end;
  end;
else
  dummy = zeros(Tmax,1);
  for t=1:Tmax,
     dummy(t)= 2*norm(R(:,:,Tmax+1+t),'fro')^2;
  end;
  dummy = cumsum(dummy,'reverse'); 
  xi_2 = zeros(T,1);
  if T<=Tmax,
     xi_2 = dummy(1:T);
  else; 
     xi_2(1:Tmax) = dummy;
  end;        
end;

%------------------------------------------------------------------------------
%  determine estimation error
%------------------------------------------------------------------------------
if T > Tmax,
   L = T-Tmax;
   dummy = zeros(M,M,2*T+1);
   dummy(:,:,L+1:L+2*Tmax+1)=R;
   R = dummy; Tmax = T; 
end;
VarE = zeros(M,M,2*T+1);
for m = 1:M,
  r_xx = squeeze(R(m,m,:)).'; 
  for mu = 1:M,
    r_yy = squeeze(R(mu,mu,:)).';
    r_xx_yy_ext = [zeros(1,T) r_xx.*conj(r_yy) zeros(1,T)];   
    Var_rh_xy2 = zeros(1,2*T+1);   
    rb_xy = squeeze(R(m,mu,:)).';
    if ComplexMode==0,
      rb_xy = zeros(size(rb_xy));
    end;
    rb_xy_ext = [zeros(1,T) rb_xy zeros(1,T)];      
    for tau = -T:T,
      W = [1:(T+1), T:-1:1] + (N-abs(tau)-T-1);    % weighting function
      W1 = r_xx_yy_ext((-T:T)+Tmax+T+1);
      W2 = rb_xy_ext((-T:T)+Tmax+T+1+tau);
      Var_rh_xy2(tau+T+1) = (W1 + W2.*fliplr(conj(W2)) )*W'/(N-abs(tau)).^2;
    end;
    VarE(m,mu,:) =  Var_rh_xy2;
  end;
end;
if VectorMode==0,
  xi_1 = sum(sum(sum(VarE)));
else
  xi_1 = zeros(T,1);
  xi_1(1) = sum(sum(VarE(:,:,T+1))); 
  xi_1(2:T) = 2*cumsum(squeeze(sum(sum(VarE(:,:,T+2:2*T+1)))));  
end;  
