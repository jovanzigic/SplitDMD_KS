function [w,e,b,niter,err,alphas] = optdmd_KS(X,t,r,varargin)
%OPTDMD Wrapper of VARPRO2 for computing the optimized DMD of data
%
%   [w,e,b] = optdmd_KS(X,t,r)
%
% Input:
%
% X - data matrix, length(t) columns with each a snapshot
%   of some time-varying system
% t - times corresponding to each snapshot
% r - rank of fit, i.e. number of exponentials to fit

% Output:
%
% w - each column is a DMD mode
% e - each entry e(i) is an eigenvalue corresponding to w(:,i)
% b - the best fit coefficient of each DMD mode

% X should be approximated by
%
% X ~ w*diag(b)*exp(e*t')
%
% if t is a column vector.
%
% [w,e,b] = optdmd_KS(X,t,r);

% Copyright Travis Askham 2017
% MIT License

%% Finite difference
[u,~,~] = svd(X,'econ');

ux1 = u'*X;
ux2 = ux1(:,2:end);
ux1 = ux1(:,1:end-1);

t1 = t(1:end-1);
t2 = t(2:end);

dx = (ux2-ux1)*diag(1./(t2-t1));
xin = (ux1+ux2)/2;

%% standard DMD
[u1,s1,v1] = svd(xin,'econ');

u1 = u1(:,1:r);
v1 = v1(:,1:r);
s1 = s1(1:r,1:r);

atilde = u1'*dx*v1/s1;

alpha_init = eig(atilde);

clear ux1 ux2 atilde t1 t2 dx xin

%% fit to all of data

m = length(t);
[Nx,~] = size(X);

% variable projection method (Levenberg-Marquardt algorithm)
[w,e,niter,err,alphas] = ...
    varpro2_KS(transpose(X),t,@varpro2expfun,@varpro2dexpfun,m,Nx,r,alpha_init);

w = w.';

%% normalize

b = sqrt(sum(abs(w).^2,1))';
inds_small = abs(b) < 10*eps(1)*max(b);
b( inds_small ) = 1.0;    
w = w*diag(1./b);
w(:,inds_small) = 0.0;
b( inds_small ) = 0.0;




