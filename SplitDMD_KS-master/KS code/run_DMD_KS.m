function [relerrDMD_r,relerr1_r,ROM_error_DMD,ROM_error_OD,tOD,...
    tDMD,niter_LM,err_LM,alphas_LM,realXdmd,errorDMD,realXod,errorOD,w1,e1,b1] = run_DMD_KS(X,r,dt,M,t,interval)
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

tic
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
tDMD = toc;

relerrDMD_r = norm(Xdmd-X1,'fro')/norm(X1,'fro');

%% Step 2b, 3b, 4b: Compute Optimized DMD modes, amplitudes and reconstruct X

% fit to unprojected data
tic
[w1,e1,b1,niter_LM,err_LM,alphas_LM] = optdmd_KS(X,t,r);

% reconstructed values
Xd1 = w1*diag(b1)*exp(e1*t);
relerr1_r = norm(Xd1-X,'fro')/norm(X,'fro');
tOD = toc;

%% Step 5: Plot DMD eigenvalues
% figure
% hold on
% scatter(real(omega),imag(omega),'bd','LineWidth',2)
% scatter(real(e1),imag(e1),'rx','LineWidth',2)
% box on
% grid on
% legend('Standard DMD eigenvalues','Optimized DMD eigenvalues')
% title(['Modes for Rank ' num2str(r), ' Approximation'])

%% Step 6: Plot Modes and Error
Nx = size(Xdmd,1);
Nt = size(Xdmd,2);

x = linspace(0,1,Nx);
t = linspace(0,interval,Nt);

% Standard DMD
[ROM_error_DMD,realXdmd,errorDMD] = simDMD...
    (Xdmd,X1,x,t,Nt,M,'Standard DMD Model',r,'z_{DMD}', 'z-z_{DMD}',interval);

% Optimized DMD
[ROM_error_OD,realXod,errorOD] = simDMD...
    (Xd1(:, 1:end-1),X1,x,t,Nt,M,'Optimized DMD Model',r,'z_{OD}', 'z-z_{OD}',interval);

end

% Simulate ROM
function [ROM_error_DMD,realXdmd,error] = simDMD(Xdmd,X1,x,t,Nt,M,name,r,model1,model2,interval)

realXdmd = real(Xdmd);
% figure
% mesh(x,t,realXdmd'); 
% view([.7 -.8 .6])
% title([name, ', order = ' , num2str(r) ] )
% xlabel('x'); ylabel('t'); zlabel(model1)
% xlim([0 1]) 
% ylim([0 interval])
% zlim([-20 20])

error = X1 - realXdmd;
% figure
% mesh(x,t,error')
% view([.7 -.8 .6])
% xlabel('x'); ylabel('t'); zlabel(model2)
% title(['Error in ', name, ', order = ' , num2str(r)] )
% xlim([0 1]) 
% ylim([0 interval])
% zlim([-20 20])

delta_t = (t(2)-t(1))/Nt;

ROM_error_DMD = 0.5*error(:,1)'*M*error(:,1);
  for i=2:Nt-1
    ROM_error_DMD = ROM_error_DMD + error(:,i)'*M*error(:,i);
  end
  ROM_error_DMD = ROM_error_DMD + 0.5*error(:,Nt)'*M*error(:,Nt);
  
  ROM_error_DMD = sqrt(delta_t*ROM_error_DMD); % (\int_0^T \| error(t,.) \|^2 dt)^{1/2}

end
