clear, close all

% KS (Jeff's code)
load kuramoto_1db_snap 
z = w_save;
t = time;

% ROM size
r_dim = 20;
dt = 1/length(t(1:end-1));

% Simulate full model
figure
mesh(x,t(1:end-1),z')
view([.7 -.8 .6])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution')
xlim([0 1]) 
ylim([0 12])
zlim([-15 15])

% Mass matrix (Jeff's code)
[M] = compute_mass_matrix(x,e_conn);

% POD
[POD] = run_POD(z,r_dim,M);

% Use the POD modes to compute a DMD bases of dimension r_dim
[relerrDMD_r,relerr1_r,ROM_error_DMD,ROM_error_optDMD,tOD,tDMD] = run_DMD(z,r_dim,dt,M,t(1:end-1));

% % Simulate POD-Galerkin reduced model
% [ROM_error_POD,relerrPOD_r] = test_rom_KS(POD,r_dim,epsilon,q1,q2,M,z,r_dim);



