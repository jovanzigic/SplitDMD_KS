clear, close all

interval = 400; % adjust parameter
re = 402.3; % adjust parameter
r_dim = 13; % adjust parameter

split = 1;

if split == 1
    fileload = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split200');
else
    fileload = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
end
load(fileload)

if split == 1
    % Get the DMD modes for each phase
    [Phi_1,omega_1,lambda_1,b_1,time_dynamics_1] = run_std_DMD(z_p1,r_dim,dt);
    [Phi_2,omega_2,lambda_2,b_2,time_dynamics_2] = run_std_DMD(z_p2,r_dim,dt);
%     [Phi_3,omega_3,lambda_3,b_3,time_dynamics_3] = run_std_DMD(z_p3,r_dim,dt);
    eps = 1e-7;
    for i = 1:(r_dim-1)
        if abs(imag(Phi_2(1,i))) < eps
            w2 = [w2(:,1:(i-1)) w2(:,(i+1):end)];
            Phi_2 = [Phi_2(:,1:(i-1)) Phi_2(:,(i+1):end)];
        end
    end
    if abs(imag(Phi_2(1,end))) < eps
        w2 = w2(:,1:(end-1));
        Phi_2 = Phi_2(:,1:(end-1));
    end
    cd_2norm = norm(w2)/norm(Phi_2);
else
    [Phi,omega,lambda,b,time_dynamics] = run_std_DMD(z,r_dim,dt);
    for i = 1:(r_dim-1)
        if abs(imag(Phi(1,i))) < eps
%             w = [w(:,1:(i-1)) w(:,(i+1):end)];
            Phi = [Phi(:,1:(i-1)) Phi(:,(i+1):end)];
        end
    end
    if abs(imag(w(1,end))) < eps
        w = w(:,1:(end-1));
        Phi = Phi(:,1:(end-1));
    end
    cd_2norm = norm(w)/norm(Phi);
end

% Save data
save(fileload)