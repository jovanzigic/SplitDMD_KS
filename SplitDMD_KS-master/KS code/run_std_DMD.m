function [Phi,omega,lambda,b,time_dynamics] = run_std_DMD(X,r,dt)
%%
% The DMD method provides spatiotemporal decomposition of data into a set
% of dynamic modes that are derived from snapshots of a system in time.
%
% INPUTS: 
% snapshots = data matrix
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X')
% M = FEM mass matrix
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstructed by Phi, omega, b

%% Step 1: 
% 'm' state snapshots from k = 1:m are arranged into two large matrices:
% X1 = snapshots from k = 1:m-1
% X2 = snapshots from k = 2:m

% Columns of X1 and X2 are state snapshots 
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

%% Step 2a: Compute the DMD modes
% the DMD modes are the eigenvectors of A

[U, S, V] = svd(X1, 'econ'); % SVD of X1
r = min(r, size(U,2)); % choose rank of reduced model

% truncate to rank-r
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

% A_r is r x r projection of A onto POD modes
A_r = U_r' * X2 * V_r / S_r; % low-rank dynamics
[W_r, D] = eig(A_r); % eigendecomposition of A_r
Phi = X2 * V_r / S_r * W_r; % DMD modes

% DMD spectra
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

%% Step 3: Compute DMD mode amplitudes b

x1 = X1(:, 1); % initial state snapshot
b = Phi\x1; % x1 = Phi*b, b = the initial amplitudes of each mode

%% Step 4: DMD reconstruction
% Xdmd = x(t), the approximate solution of future states

time_dynamics = zeros(r, length(X1));
t_vec = (0:length(X1)-1)*dt; % time vector
for i = 1:length(X1)
    time_dynamics(:,i) = (b.*exp(omega*t_vec(i)));
end
Xdmd = Phi * time_dynamics;

