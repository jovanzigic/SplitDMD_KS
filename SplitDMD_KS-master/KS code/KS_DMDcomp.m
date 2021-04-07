% Chaotic Case Metric for Pattern Reconstruction
% || OD modes ||_2 / || DMD modes ||_2 for each split

clear, close all

interval = 400; % adjust parameter
re = 402.3; % adjust parameter
r_dim = 13; % adjust parameter

split = 1;

if split == 1
    fileload1 = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split4_error140770');
else
    n_split = 1;
    fileload1 = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
end
load(fileload1)

% Obtain modes for each split
Phi = cell(1,n_split);
cd_2norm = zeros(1,n_split);
w1 = w;
for k = 1:n_split
    [Phi{1,k},~,~,~,~] = run_std_DMD(z(:,p(k):p(k+1)),r_dim,dt);
    Phimat = cell2mat(Phi(1,k));
    wmat = cell2mat(w1(1,k));
    ct = 0;
    % DMD Modes
    eps = 1e-15;
    Phi1 = Phimat;
    for i = 1:(r_dim-1)
        if abs(imag(Phimat(1,i))) < eps
%             wmat = [wmat(:,1:(i-1)) wmat(:,(i+1):end)];
            Phi1 = [Phi1(:,1:(i-1-ct)) Phi1(:,(i+1-ct):end)];
            ct = ct + 1;
        end
    end
    Phimat = Phi1;
    if abs(imag(Phimat(1,end))) < eps
%         wmat = wmat(:,1:(end-1));
        Phimat = Phimat(:,1:(end-1));
    end
    ct = 0;
    % OD Modes
    eps = 1e-15;
    w2 = wmat;
    for i = 1:(r_dim-1)
        if abs(imag(wmat(1,i))) < eps
            w2 = [w2(:,1:(i-1-ct)) w2(:,(i+1-ct):end)];
            ct = ct + 1;
        end
    end
    wmat = w2;
    if (abs(imag(wmat(1,end))) < eps) || (size(Phimat,2)+1 == size(wmat,2))
        wmat = wmat(:,1:(end-1));
    end
    Phi{1,k} = Phimat;
    w1{1,k} = wmat;
    if size(Phimat,2) == size(wmat,2)
        cd_2norm(k) = norm(wmat)/norm(Phimat);
    else
        cd_2norm(k) = 1;
    end
end

% Save data
save(fileload1)