clear, close all

% choose rank loop
N = 10;
space = 1;
n = N/space;

% choose # of reynolds loops
ct = 1;

% array of rank approximations for a reynolds number
dim = zeros(1,n);
for i = 1:n
    dim(1,i) = space*i + 4 ;
end

% collection of errors for plotting
% PODerrF_full = zeros(ct,n);
DMDerrF_full = zeros(ct,n);
ODerrF_full = zeros(ct,n);
% % PODerr2_full = zeros(ct,n);
% DMDerr2_full = zeros(ct,n);
% ODerr2_full = zeros(ct,n);

iter1 = 0;
timeOD = zeros(1,ct*n);
timeDMD = zeros(1,ct*n);
% timePOD = zeros(1,ct*n);
timespace = linspace(1,ct*n,ct*n);

% KS (Jeff's code)
re = 13.2;
load kuramoto_1db_snap_L1380
z = w_save;
t = time;

% initialize array of errors
% PODerrF = zeros(1,n);
DMDerrF = zeros(1,n);
ODerrF = zeros(1,n);

% For multiple Re loops
iter1 = iter1 + 1;

for rk = 1:n

    % ROM size (the number of basis functions to use)
    r_dim  = dim(1,rk);
    dt = 1/length(t(1:end-1));

%     % Simulate full model
%     figure
%     mesh(x,t(1:end-1),z')
%     view([.7 -.8 .6])
%     pause(0.001)
%     xlabel('x'); ylabel('t'); zlabel('z')
%     title('Finite Element Solution')
%     xlim([0 1]) 
%     ylim([0 12])
%     zlim([-15 15])

    % Mass matrix (Jeff's code)
    [M] = compute_mass_matrix(x,e_conn);

    % POD
    [POD] = run_POD(z,r_dim,M);

    % Use the POD modes to compute a DMD bases of dimension r_dim
    [relerrDMD_r,relerr1_r,~,~,tOD,tDMD,niter_LM,err_LM,alphas_LM] ...
        = run_DMD_KS(z,r_dim,dt,M,t(1:end-1));

    % % Simulate POD-Galerkin reduced model
    % [ROM_error_POD,relerrPOD_r] = test_rom_KS(POD,r_dim,epsilon,q1,q2,M,z,r_dim);

    % Collect error norms
    % PODerrF(1,rk) = relerrPOD_r;
    DMDerrF(1,rk) = relerrDMD_r;
    ODerrF(1,rk) = relerr1_r;
    
    % For 1 Re, 6 rk
    timeOD(1, rk) = tOD;
    timeDMD(1, rk) = tDMD;
    % timePOD(1, (iter2 - 1)*5 + iter1) = tPOD;
    
end

% PODerrF_full(iter1,:) = PODerrF;
DMDerrF_full(iter1,:) = DMDerrF;
ODerrF_full(iter1,:) = ODerrF;

% Error Plots for single Re
figure('Name','Relative Error of ROM')
semilogy(dim,DMDerrF,'b-+',dim,ODerrF,'r-o')
title(['Relative ROM Error $ = \frac{\|Z-Z_r \|_F}{\|Z \|_F} $, Re = ', num2str(re)], 'Interpreter','latex')
legend('DMD', 'OD')
xlim([dim(1,1) dim(1,n)])
xlabel('Dimension (Order of Approximation)')
ylabel('Size of Error')
ylim([0 0.1])

% Time plot
figure('Name','Computation Time')
plot(dim,timeDMD(1,:),'b-+',dim,timeOD(1,:),'m-x')
title('Computation Time')
legend('DMD', 'OD', 'Location','northwest')
xlim([dim(1,1) dim(1,n)])
xlabel('Dimension (Order of Approximation)')
ylabel('Seconds')
