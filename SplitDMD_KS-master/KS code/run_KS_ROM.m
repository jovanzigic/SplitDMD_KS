clear, close all


%% Turn Splitting On or Off
split = 1; % adjust parameter (1 on, 0 off)
sens = 1;

%% Choose Rank of Approximation

% % N = space implies only one test
% N = 6; % adjust parameter
% space = 2; % adjust parameter
N = 13; space = 13; % adjust parameter
n = N/space;

% array of rank approximations
dim = zeros(1,n);
for i = 1:n
    dim(1,i) = space*i  ;
end

%% Initialize

% choose # of reynolds loops
ct = 1;

% collection of errors for plotting
DMDerrF_full = zeros(ct,n);
ODerrF_full = zeros(ct,n);

% collection of computation time
iter1 = 0;
timeOD = zeros(1,ct*n);
timeDMD = zeros(1,ct*n);
timespace = linspace(1,ct*n,ct*n);

% initialize array of errors 
DMDerrF = zeros(1,n);
ODerrF = zeros(1,n);

% For multiple Re loops
iter1 = iter1 + 1;

%% Load Data (double check 'dim' above)

% KS (Jeff's code)
interval = 400; % adjust parameter, time interval length
re = 402.35; % adjust parameter, bifurcation parameter

if sens == 1
    fileload = strcat('kuramoto_1db_snap_L',int2str(100*re),'_',int2str(interval),'_sens1');
else
    fileload = strcat('kuramoto_1db_snap_L',int2str(100*re),'_',int2str(interval));
end
load(fileload)
Z = w_save;
t = time;

% Simulate full model
figure
% 4000 case has 5x larger steps
if interval == 4000
    mesh(x,t(1:end),Z')
else
    mesh(x,t(1:end-1),Z')
end
% view([.7 -.8 .6])
view([0 0 1])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution')
xlim([0 1]) 
ylim([0 interval])

%% Split the Data
if split == 1
    s = 10; % adjust parameter, number of initial splits
    t_spliti = zeros(1,s);
    iter_spliti = 1;
    ni = 10;
    t_split = zeros(1,1);
    while size(t_split,2) < s
        tlist = zeros(1,1);
        for i = 1:ni
            [~,t_split,iter_split] = run_nsplit(Z,Z,s,t_spliti,iter_spliti);
            for k = 1:size(t_split,2)
                tlist(end+1) = t_split(k);
            end
        end
        t_split = zeros(1,1);
        init = true;
        for i = 1:size(tlist,2)
            if nnz(tlist==tlist(i)) > (0.6*ni)
                if init
                    init = false;
                    t_split(1) = tlist(i);
                else
                    t_split(end+1) = tlist(i);
                end
            end
        end
        t_split = sort(t_split);
        t_split = unique(t_split.','rows').';
        n_split = size(t_split,2); 
    end

%     % shift split lines
%     l1 = length(t_split)-1;
%     
%     for k=1:l1
%         t_split(k) = t_split(k) + 150;
%     end

    % splits
    p = ones(1,n_split+1);
    for k = 1:n_split
        p(1,k+1) = t_split(k);
    end

    % Simulate split k model
    for k = 1:n_split
        figure
        mesh(x,t(p(k):p(k+1)),Z(:,p(k):p(k+1))')
        % view([.7 -.8 .6])
        view([0 0 1])
        pause(0.001)
        xlabel('x'); ylabel('t'); zlabel('z')
        title(['Finite Element Solution, Split ', num2str(k)], 'Interpreter','latex')
        xlim([0 1])
        ylim([0 interval])
    end
end

%% Future State Prediction (Incomplete)

% % Extrapolating Data
% interval2 = 4000; % adjust parameter
% % load data
% fileload2 = strcat('kuramoto_1db_snap_L',int2str(100*re),'_',int2str(interval2));
% load(fileload2)
% z2 = w_save;
% t2 = time;

% % for predicting t=4000
% % p1_2 = p1/5;
% p1_2 = p1/5;
% z2_p2 = z2(:,p1_2:length(z2));

% % Simulate split 2 model t=4000
% figure
% mesh(x,t2(p1_2:end),z2_p2')
% % view([.7 -.8 .6])
% view([0 0 1])
% pause(0.001)
% xlabel('x'); ylabel('t'); zlabel('z')
% title('Finite Element Solution, Phase 2')
% xlim([0 1]) 
% ylim([0 interval2])
% % zlim([-20 20])

% % Simulate full model t=4000
% figure
% % 4000 case has 5x larger steps
% if interval2 == 4000
%     mesh(x,t2(1:end),z2')
% else
%     mesh(x,t2(1:end-1),z2')
% end
% % view([.7 -.8 .6])
% view([0 0 1])
% pause(0.001)
% xlabel('x'); ylabel('t'); zlabel('z')
% title('Finite Element Solution')
% xlim([0 1]) 
% ylim([0 interval2])
% % zlim([-20 20])

%% Run DMD
for rk = 1:n
    
      if split == 1
          
        % Mass matrix (Jeff's code)
        [M] = compute_mass_matrix(x,e_conn);
        
        % store data
        realXdmd = cell(1,n_split);
        errorDMD = cell(1,n_split);
        realXod = cell(1,n_split);
        errorOD = cell(1,n_split);
        
        relerrDMD_r = zeros(1,n_split);
        relerr1_r = zeros(1,n_split);
        tOD = zeros(1,n_split);
        tDMD = zeros(1,n_split);
        niter_LM = zeros(1,n_split);
        
        err_LM = cell(1,n_split);
        alphas_LM = cell(1,n_split);
        w = cell(1,n_split);
        e = cell(1,n_split);
        b = cell(1,n_split);
        
        for k = 1:n_split
          
            % phase k

            % ROM size (the number of basis functions to use)
            r_dim  = dim(1,rk);
            dt = 1/length(t(p(k):p(k+1)));

            % Compute a DMD bases of dimension r_dim
            [relerrDMD_r(1,k),relerr1_r(1,k),~,~,...
                tOD(1,k),tDMD(1,k),niter_LM(1,k),...
                err_LM{1,k},alphas_LM{1,k},...
                realXdmd{1,k},errorDMD{1,k},realXod{1,k},errorOD{1,k},...
                w{1,k},e{1,k},b{1,k} ] ...
                = run_DMD_KS(Z(:,p(k):p(k+1)),r_dim,dt,M,t(p(k):p(k+1)),interval);

        end
        
        % convert data
        realXdmd = cell2mat(realXdmd);
        errorDMD = cell2mat(errorDMD);
        realXod = cell2mat(realXod);
        errorOD = cell2mat(errorOD);
        
        relerrDMD_r = sum(relerrDMD_r);
        relerr1_r = sum(relerr1_r);
        tOD = sum(tOD);
        tDMD = sum(tDMD);
        
      else
            
        % ROM size (the number of basis functions to use)
        r_dim  = dim(1,rk);
        dt = 1/length(t(1:end-1));

        tic
        % Mass matrix (Jeff's code)
        [M] = compute_mass_matrix(x,e_conn);

        % Compute a DMD bases of dimension r_dim
        [relerrDMD_r,relerr1_r,~,~,tOD,tDMD,niter_LM,err_LM,alphas_LM,...
            realXdmd,errorDMD,realXod,errorOD,w,e,b] ...
            = run_DMD_KS(Z,r_dim,dt,M,t(1:end-1),interval);
        
      end

    % Simulate DMD ROM and collect L2 error
    [errorL2_DMD,errorL2_OD] = sim_DMD_KS(realXdmd,errorDMD,realXod,errorOD,r_dim,interval,M);
    
%     % Extrapolation
%     Xd10 = w*diag(b)*exp(e*t2);
%     relerr10_r = norm(Xd10-z2,'fro')/norm(z2,'fro');

%     % Simulate Extrapolation
%     realXd10 = real(Xd10);
%     figure
%     mesh(x,t2(p1_2:end),realXd10'); 
%     view([0 0 1])
%     title(['Optimized DMD Model, order = ' , num2str(r_dim) ] )
%     xlabel('x'); ylabel('t'); zlabel('z_{OD}')
%     xlim([0 1]) 
%     ylim([0 interval2])
%     
%     % Simulate Extrapolation Error
%     errorXd10 = z2_p2 - realXd10;
%     figure
%     mesh(x,t2(p1_2:end),errorXd10')
%     %view([.7 -.8 .6])
%     view([0 0 1])
%     xlabel('x'); ylabel('t'); zlabel('z-z_{OD}')
%     title(['Error in Optimized DMD Model, order = ' , num2str(r_dim)] )
%     xlim([0 1]) 
%     ylim([0 interval2])
    
%     % Compute L2 error of Extrapolation
%     ROM_error_Xd10 = 0.5*errorXd10(:,1)'*M*errorXd10(:,1);
%     for i=2:Nt-1
%     ROM_error_Xd10 = ROM_error_Xd10 + errorXd10(:,i)'*M*errorXd10(:,i);
%     end
%     ROM_error_Xd10 = ROM_error_Xd10 + 0.5*errorXd10(:,Nt)'*M*errorXd10(:,Nt);
% 
%     ROM_error_Xd10 = sqrt(delta_t*ROM_error_Xd10); 
    
    % Collect error norms
    DMDerrF(1,rk) = relerrDMD_r;
    ODerrF(1,rk) = relerr1_r;
    
    % Time
    timeOD(1, rk) = tOD;
    timeDMD(1, rk) = tDMD;
    % timePOD(1, (iter2 - 1)*5 + iter1) = tPOD;

    % Save data
    if split > 0
        filesave = strcat('ROM_KS_data_t',...
            int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),...
            '_split',int2str(n_split),'_error',int2str(10000*errorL2_OD));
    else
        filesave = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
    end
    save(filesave)
    
end

%% Plot Error and Time Results for Multiple Tests

if n > 1

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
