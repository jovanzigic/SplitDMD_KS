clear, close all

% current change : rank loop, dmd split, name = split2
split = 1; % adjust parameter

% choose rank loop
% N = 6;
% space = 2;
N = 13; space = 13;
n = N/space;

% choose # of reynolds loops
ct = 1;

% array of rank approximations
dim = zeros(1,n);
for i = 1:n
    dim(1,i) = space*i  ;
end

% collection of errors for plotting
% PODerrF_full = zeros(ct,n);
DMDerrF_full = zeros(ct,n);
ODerrF_full = zeros(ct,n);
% % PODerr2_full = zeros(ct,n);
% DMDerr2_full = zeros(ct,n);
% ODerr2_full = zeros(ct,n);

% collection of computation time
iter1 = 0;
timeOD = zeros(1,ct*n);
timeDMD = zeros(1,ct*n);
% timePOD = zeros(1,ct*n);
timespace = linspace(1,ct*n,ct*n);

% KS (Jeff's code)
interval = 120; % adjust parameter
re = 402.3; % adjust parameter
% load data (double check 'dim' above)
fileload = strcat('kuramoto_1db_snap_L',int2str(100*re),'_',int2str(interval));
load(fileload)
z = w_save;
t = time;

% Extrapolating Data
interval2 = 4000; % adjust parameter
% load data
fileload2 = strcat('kuramoto_1db_snap_L',int2str(100*re),'_',int2str(interval2));
load(fileload2)
z2 = w_save;
t2 = time;

% initialize array of errors 
% PODerrF = zeros(1,n);
DMDerrF = zeros(1,n);
ODerrF = zeros(1,n);

% For multiple Re loops
iter1 = iter1 + 1;

% Phases
p1 = (200/400)*length(z);
% p2 = length(z)-p1;
p2 = (400/400)*length(z);
% p3 = (400/400)*length(z);
z_p1 = z(:,1:p1);
% z_p2 = z(:,p1:length(z));
z_p2 = z(:,p1:p2);
% z_p3 = z(:,p2:p3);

% p1_2 = p1/5;
p1_2 = p1/5;
z2_p2 = z2(:,p1_2:length(z2));

% Simulate phase 1 model
figure
mesh(x,t(1:p1),z_p1')
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution, Phase 1')
xlim([0 1]) 
ylim([0 interval])
% zlim([-20 20])

% Simulate phase 2 model
figure
mesh(x,t(p1:p2),z_p2')
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution, Phase 2')
xlim([0 1]) 
ylim([0 interval])
% zlim([-20 20])

% % Simulate phase 3 model
% figure
% mesh(x,t(p2:p3),z_p3')
% % view([.7 -.8 .6])
% view([0 0 1])
% pause(0.001)
% xlabel('x'); ylabel('t'); zlabel('z')
% title('Finite Element Solution, Phase 3')
% xlim([0 1]) 
% ylim([0 interval])
% % zlim([-20 20])

% Simulate phase 2 model t=4000
figure
mesh(x,t2(p1_2:end),z2_p2')
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution, Phase 2')
xlim([0 1]) 
ylim([0 interval2])
% zlim([-20 20])

% Simulate full model
figure
% 4000 case has 5x larger steps
if interval == 4000
    mesh(x,t(1:end),z')
else
    mesh(x,t(1:end-1),z')
end
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution')
xlim([0 1]) 
ylim([0 interval])
% zlim([-20 20])

% Simulate full model
figure
% 4000 case has 5x larger steps
if interval2 == 4000
    mesh(x,t2(1:end),z2')
else
    mesh(x,t2(1:end-1),z2')
end
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution')
xlim([0 1]) 
ylim([0 interval2])
% zlim([-20 20])

for rk = 1:n
    
      if split == 1
          
        % phase 1
    
        % ROM size (the number of basis functions to use)
        r_dim  = dim(1,rk);
        dt = 1/length(t(1:p1));

        % Mass matrix (Jeff's code)
        [M] = compute_mass_matrix(x,e_conn);

        % POD
        [POD1] = run_POD(z_p1,r_dim,M);

        % Use the POD modes to compute a DMD bases of dimension r_dim
        [relerrDMD_r_p1,relerr1_r_p1,~,~,tOD_p1,tDMD_p1,niter_LM_p1,err_LM_p1,alphas_LM_p1,realXdmd_p1,...
            errorDMD_p1,realXod_p1,errorOD_p1,w1,e1,b1] ...
            = run_DMD_KS(z_p1,r_dim,dt,M,t(1:p1),interval);
    
%         % phase 2
%     
%         % ROM size (the number of basis functions to use)
%         dt = 1/length(t(p1:end-1));
% 
%         % POD
%         [POD2] = run_POD(z_p2,r_dim,M);
% 
%         % Use the POD modes to compute a DMD bases of dimension r_dim
%         [relerrDMD_r_p2,relerr1_r_p2,~,~,tOD_p2,tDMD_p2,niter_LM_p2,err_LM_p2,alphas_LM_p2,realXdmd_p2,...
%             errorDMD_p2,realXod_p2,errorOD_p2,w2,e2,b2] ...
%             = run_DMD_KS(z_p2,r_dim,dt,M,t(p1:end-1),interval);

        % phase 2
    
        % ROM size (the number of basis functions to use)
        dt = 1/length(t(p1:p2));

        % POD
        [POD2] = run_POD(z_p2,r_dim,M);

        % Use the POD modes to compute a DMD bases of dimension r_dim
        [relerrDMD_r_p2,relerr1_r_p2,~,~,tOD_p2,tDMD_p2,niter_LM_p2,err_LM_p2,alphas_LM_p2,realXdmd_p2,...
            errorDMD_p2,realXod_p2,errorOD_p2,w2,e2,b2] ...
            = run_DMD_KS(z_p2,r_dim,dt,M,t(p1:p2),interval);
        
%         % phase 3
%     
%         % ROM size (the number of basis functions to use)
%         dt = 1/length(t(p2:p3));
% 
%         % POD
%         [POD3] = run_POD(z_p3,r_dim,M);
% 
%         % Use the POD modes to compute a DMD bases of dimension r_dim
%         [relerrDMD_r_p3,relerr1_r_p3,~,~,tOD_p3,tDMD_p3,niter_LM_p3,err_LM_p3,alphas_LM_p3,realXdmd_p3,...
%             errorDMD_p3,realXod_p3,errorOD_p3,w3,e3,b3] ...
%             = run_DMD_KS(z_p3,r_dim,dt,M,t(p2:p3),interval);
%                 
        % summation for 2 phases
        realXdmd = [realXdmd_p1 realXdmd_p2];
        errorDMD = [errorDMD_p1 errorDMD_p2];
        realXod = [realXod_p1 realXod_p2];
        errorOD = [errorOD_p1 errorOD_p2];
        relerrDMD_r = relerrDMD_r_p1 + relerrDMD_r_p2;
        relerr1_r = relerr1_r_p1 + relerr1_r_p2;
        tOD = tOD_p1 + tOD_p2;
        tDMD = tDMD_p1 + tDMD_p2;

%         % summation for 3 phases
%         realXdmd = [realXdmd_p1 realXdmd_p2 realXdmd_p3];
%         errorDMD = [errorDMD_p1 errorDMD_p2 errorDMD_p3];
%         realXod = [realXod_p1 realXod_p2 realXod_p3];
%         errorOD = [errorOD_p1 errorOD_p2 errorOD_p3];
%         relerrDMD_r = relerrDMD_r_p1 + relerrDMD_r_p2+ relerrDMD_r_p3;
%         relerr1_r = relerr1_r_p1 + relerr1_r_p2 + relerr1_r_p3;
%         tOD = tOD_p1 + tOD_p2 + tOD_p3;
%         tDMD = tDMD_p1 + tDMD_p2 + tDMD_p3;

        
        % Extrapolation
        Xd10 = w2*diag(b2)*exp(e2*t2(p1_2:end));
        relerr10_r = norm(Xd10-z2_p2,'fro')/norm(z2_p2,'fro');
        
      else
            
        % ROM size (the number of basis functions to use)
        r_dim  = dim(1,rk);
        dt = 1/length(t(1:end-1));

        tic
        % Mass matrix (Jeff's code)
        [M] = compute_mass_matrix(x,e_conn);
        tM = toc;

        % POD
        [POD] = run_POD(z,r_dim,M);

        % Use the POD modes to compute a DMD bases of dimension r_dim
        [relerrDMD_r,relerr1_r,~,~,tOD,tDMD,niter_LM,err_LM,alphas_LM,realXdmd,errorDMD,realXod,errorOD,w,e,b] ...
            = run_DMD_KS(z,r_dim,dt,M,t(1:end-1),interval);
    
        % Extrapolation
        Xd10 = w*diag(b)*exp(e*t2);
        relerr10_r = norm(Xd10-z2,'fro')/norm(z2,'fro');
        
      end

    % Simulate DMD ROM
    [errorL2_DMD,errorL2_OD] = sim_DMD_KS(realXdmd,errorDMD,realXod,errorOD,r_dim,interval,M);
    
    % Simulate Extrapolation
    realXd10 = real(Xd10);
    figure
    mesh(x,t2(p1_2:end),realXd10'); 
    view([0 0 1])
    title(['Optimized DMD Model, order = ' , num2str(r_dim) ] )
    xlabel('x'); ylabel('t'); zlabel('z_{OD}')
    xlim([0 1]) 
    ylim([0 interval2])
    
    % Simulate Extrapolation Error
    errorXd10 = z2_p2 - realXd10;
    figure
    mesh(x,t2(p1_2:end),errorXd10')
    %view([.7 -.8 .6])
    view([0 0 1])
    xlabel('x'); ylabel('t'); zlabel('z-z_{OD}')
    title(['Error in Optimized DMD Model, order = ' , num2str(r_dim)] )
    xlim([0 1]) 
    ylim([0 interval2])
    
%     % Compute L2 error of Extrapolation
%     ROM_error_Xd10 = 0.5*errorXd10(:,1)'*M*errorXd10(:,1);
%     for i=2:Nt-1
%     ROM_error_Xd10 = ROM_error_Xd10 + errorXd10(:,i)'*M*errorXd10(:,i);
%     end
%     ROM_error_Xd10 = ROM_error_Xd10 + 0.5*errorXd10(:,Nt)'*M*errorXd10(:,Nt);
% 
%     ROM_error_Xd10 = sqrt(delta_t*ROM_error_Xd10); 
    
    % % Simulate POD-Galerkin reduced model
    % [ROM_error_POD,relerrPOD_r] = test_rom_KS(POD,r_dim,epsilon,q1,q2,M,z,r_dim);

    % Collect error norms
    % PODerrF(1,rk) = relerrPOD_r;
    DMDerrF(1,rk) = relerrDMD_r;
    ODerrF(1,rk) = relerr1_r;
    
    % Time
    timeOD(1, rk) = tOD;
    timeDMD(1, rk) = tDMD;
    % timePOD(1, (iter2 - 1)*5 + iter1) = tPOD;

    % Save data
    if split > 0
        filesave = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split200');
    else
        filesave = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
    end
    save(filesave)
    
end

if n > 1

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
    plot(dim,timeDMD(1,:),'b-+',dim,timeOD(1,:),'r-o')
    title('Computation Time')
    legend('DMD', 'OD', 'Location','northwest')
    xlim([dim(1,1) dim(1,n)])
    xlabel('Dimension (Order of Approximation)')
    ylabel('Seconds')

end
