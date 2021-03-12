function [b,alpha,niter,err,alphas] = varpro2_KS(y,t,phi,dphi, ...
    m,is,ia,alpha_init)
%VARPRO2 Variable projection algorithm for multivariate data
%
% Attempts a fit of the columns of y as linear combinations
% of the columns of phi(alpha,t), i.e.
%
% y_k = sum_j=1^n b_jk phi_j(alpha,t)
%
% Note that phi(alpha,t) is a matrix of dimension
% m x n where m is length (t) and n is number of columns.
%
% phi_j(alpha,t) is the jth column
% y_k is the kth column of the data
%
% Input:
%
% y - M x IS matrix of data
% t - M vector of sample times
% phi(alpha,t) - M x N matrix (or sparse matrix) valued 
%              function with input alpha
% dphi(alpha,t,i) - M x N matrix (or sparse matrix) valued
%                 function of alpha: returns the derivative 
%                 of the entries of phi with respect to the 
%                 ith component of alpha
% m - integer, number of rows of data/number of sample times
% is - integer, number of columns of data .. number of 
%      functions to fit
% ia - integer, dimension of alpha
% alpha_init - initial guess for alpha

% Optimization parameters:
% 'lambda0' (1.0) --- lambda0 is the initial 
%   value used for the regularization parameter
%   lambda in the Levenberg method (a larger
%   lambda makes it more like gradient descent)
%
% 'maxlam' (52) --- maxlam is the maximum number 
%   of steps used in the inner Levenberg loop,
%   i.e. the number of times you increase lambda
%   before quitting
%
% 'lamup' (2.0) --- lamup is the factor by which
%   you increase lambda when searching for an 
%   appropriate step
%
% 'maxiter' (50) --- the maximum number of outer
%   loop iterations to use before quitting

% 'tol' (1.0e-6) --- the tolerance for the relative
%   error in the residual, i.e. the program will
%   terminate if 
%       norm(y-Phi(alpha)*b,'fro')/norm(y,'fro') < tol
%   is achieved.
%
% 'eps_stall' (1.0e-12) --- the tolerance for detecting 
%   a stall. If err(iter-1)-err(iter) < eps_stall*err(iter-1)
%   then a stall is detected and the program halts.

% Output:
%
% b - N x IS matrix of coefficients .. each column gives
%     the coefficients for one of the functions (columns
%     of data) corresponding to the best fit
% alpha - N vector of values of alpha for best fit
% niter - number of iterations of the Marquardt algorithm
% err - the error for each iteration of the algorithm

% Copyright (c) 2018 Travis Askham
% Available under the MIT license
%
% References: 
% - Extensions and Uses of the Variable Projection 
% Algorithm for Solving Nonlinear Least Squares Problems by 
% G. H. Golub and R. J. LeVeque ARO Report 79-3, Proceedings 
% of the 1979 Army Numerical Analsysis and Computers Conference
% - "Variable projection for nonlinear least squares problems." 
% Computational Optimization and Applications 54.3 (2013): 579-593. 
% by Dianne P. O’Leary and Bert W. Rust. 

% set optimization parameters

lambda0 = 1;
maxlam = 52;
lamup = 2;
maxiter = 50;
tol = 1e-6;
eps_stall = 1e-12;

% initialize values

alpha = alpha_init;
alphas = zeros(length(alpha),maxiter);
djacmat = zeros(m*is,ia);
rhstemp = zeros(m*is,1);
err = zeros(maxiter,1);
scales = zeros(ia,1);

rjac = zeros(2*ia,ia);

phimat = phi(alpha,t);
[U,S,V] = svd(phimat,'econ');
sd = diag(S);
tolrank = m*eps;
irank = sum(sd > tolrank*sd(1));
U = U(:,1:irank);
S = S(1:irank,1:irank);
V = V(:,1:irank);
b = phimat\y;
res = y - phimat*b;
errlast = 0.5*(norm(res,'fro')^2);

for iter = 1:maxiter
  
  % build jacobian matrix, looping over alpha indices
  for j = 1:ia
    dphitemp = dphi(alpha,t,j);
    djaca = (dphitemp - sparse(U*(sparse(U'*dphitemp))))*b;
    
	% use full expression for Jacobian
    djacb = U*(S\(V'*(sparse(dphitemp'*res))));
    djacmat(1:m*is,j) = (djaca(:) + djacb(:));

	% the scales give the "marquardt" part of the algo.
    scales(j) = 1;
    scales(j) = min(norm(djacmat(1:m*is,j)),1);
    scales(j) = max(scales(j),1e-6);
  end

  % loop to determine lambda (lambda gives the "levenberg" part)
	% pre-compute components that don't depend on 
	% step-size parameter (lambda)
  
  % get pivots and lapack style qr for jacobian matrix
  rhstemp(1:m*is) = res(:);

  g = djacmat'*rhstemp;
  
  [qout,djacout,jpvt] = qr(djacmat,0);
  ijpvt = 1:ia;
  ijpvt(jpvt) = ijpvt;
  rjac(1:ia,:) = triu(djacout(1:ia,:));
  rhstop = qout'*rhstemp;
  scalespvt = scales(jpvt(1:ia)); % permute scales appropriately...
  rhs = [rhstop(1:ia); zeros(ia,1)]; % transformed right hand side
  
  % check if current step size or shrunk version works
    % get step

  rjac(ia+1:2*ia,:) = lambda0*diag(scalespvt);

  delta0 = rjac\rhs;
  delta0 = delta0(ijpvt); % unscramble solution
  
  % new alpha guess
  alpha0 = alpha + delta0;
  
  % corresponding residual
  phimat = phi(alpha0,t);
%   for i = 1:size(phimat,1)
%       for j = 1:size(phimat,2)
%           if phimat(i,j) < 10e-06
%             phimat(i,j) = 0;
%           end
%       end
%   end
  b0 = phimat\y;
  res0 = y-phimat*b0;
  err0 = 0.5*(norm(res0,'fro')^2);
  
  % new update rule: check
  % predicted improvement vs actual improvement
  act_impr = errlast-err0;
  pred_impr = real(0.5*delta0'*(g));
  impr_rat = act_impr/pred_impr;
  
  if (err0 < errlast)
      % new version
      % rescale lambda based on
      % actual vs pred improvement    
      lambda0 = lambda0*max(1.0/3.0,1-(2*impr_rat-1)^3);
      alpha = alpha0;
      errlast = err0;
      b = b0;
      res = res0;   
  else
	% if not, increase lambda until something works
	% this makes the algorithm more and more like gradient descent
    for j = 1:maxlam
      
      lambda0 = lambda0*lamup;
      rjac(ia+1:2*ia,:) = lambda0*diag(scalespvt);

      delta0 = rjac\rhs;
      delta0 = delta0(ijpvt); % unscramble solution
      
      alpha0 = alpha + delta0;

      phimat = phi(alpha0,t);
      b0 = phimat\y;
      res0 = y-phimat*b0;
      err0 = 0.5*(norm(res0,'fro')^2);
      
      if (err0 < errlast) 
        break
      end

    end
    
    if (err0 < errlast) 
      alpha = alpha0;
      errlast = err0;
      b = b0;
      res = res0;
    else
      
	  % no appropriate step length found
      niter = iter;
      err(niter) = errlast;
      return
    end
  end
  
  alphas(:,iter) = alpha;
  
  err(iter) = errlast;
  
  if (errlast < tol)
	% tolerance met    
    niter = iter;
    return;
  end
  
  if (iter > 1)
    if (err(iter-1)-err(iter) < eps_stall*err(iter-1))
	  % stall detected      
      niter = iter;
      return;
    end
  end
  
  phimat = phi(alpha,t);
  [U,S,V] = svd(phimat,'econ');
  sd = diag(S);
  irank = sum(sd > tolrank*sd(1));
  U = U(:,1:irank);
  S = S(1:irank,1:irank);
  V = V(:,1:irank);
  
end

% failed to meet tolerance in maxiter steps
niter = maxiter;

end


