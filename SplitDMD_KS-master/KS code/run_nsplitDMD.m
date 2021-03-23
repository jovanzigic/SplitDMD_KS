function [n_split,t_split,iter] = run_nsplitDMD(data,z,s,t_split,iter,r,t)
% starting from initial split in two parts, the data is split further
% if the intervals have proportionally distinct data ranges 

% z = data
% x = chosen row of data (one space over time)
% size(z,1) = spatial data
% size(z,2) = temporal data
% s (input), n_split (output) = number of splits
% t_split = temporal split location, initially = zeros(1,s)
% iter = count of iteration, initially 1
% p = column cutoff for split of z
% og = original size of time interval

% choose epsilon tolerance for z (fraction of data range)
% choose delta tolerance for x (fraction of time interval)
    
    n_split = s;

%     % choose random space across time
%     x_rand = randi(size(z,1),1,1);
% 
%     x = x_rand;

    eps = round((1/5)*range(data(:,:)));
    del = round((1/5)*size(data(:,:),2));

    % initialize 
    p = ones(1,s+1);
    stall = true;

    if t_split == zeros(1,size(t_split,2))
        for k = 1:s
            t_split(1,k) = round((k/s)*size(z,2));
        end
    end

    t_split = sort(t_split);
    t_split = unique(t_split.','rows').';

    for k = 1:s
        p(1,k+1) = round((t_split(1,k)/t_split(1,size(t_split,2)))*size(z,2));
    end

    % Get the DMD modes for each phase
    Phi = cell(1,s);
    eps = 1e-20;
    Phi_norm = zeros(1,s);
    for k = 1:s 
        dt = 1/length(t(p(k):p(k+1)));
        [Phi{1,k},~,~,~,~] = run_std_DMD(z(:,p(k):p(k+1)),r,dt);
        Phimat = cell2mat(Phi(1,k));
        ct = 0;
        Phi1 = Phimat;
        for i = 1:(r-1)
            if abs(imag(Phimat(1,i))) < eps
                Phi1 = [Phi1(:,1:(i-1-ct)) Phi1(:,(i+1-ct):end)];
                ct = ct + 1;
            end
        end
        Phimat = Phi1;
        if abs(imag(Phimat(1,end))) < eps
            Phimat = Phimat(:,1:(end-1));
        end
        Phi{1,k} = Phimat;
        Phi_norm(k) = norm(Phimat(:,1:2));
    end
    
%     g1 = (cell2mat(Phi(1,1)));
%     g2 = (cell2mat(Phi(1,2)));
%     g3 = (cell2mat(Phi(1,3)));
%     g4 = (cell2mat(Phi(1,4)));
    
%     cd_2norm = norm(cell2mat(Phi(1,2)))/norm(cell2mat(Phi(1,1)));
    cd_2norm = norm(cell2mat(Phi(1,2)))/norm(cell2mat(Phi(1,1)));
    
    s_1 = 1; 
    s_2 = 1; 
    s_k = s;    
    
    for k = 1:(s-1)
        if (abs( max(z(:,p(1,k+1):p(1,k+2))) - max(z(:,p(1,k):p(1,k+1))) ) > eps ...
                || abs( min(z(:,p(1,k+1):p(1,k+2))) - min(z(:,p(1,k):p(1,k+1))) ) > eps)
            if abs(p(1,k+1) - p(1,k)) > del
                stall = false;
                t_1 = [ t_split round((p(1,k+1) + p(1,k))/2) ];
                [s_1,t_1,iter1] = run_nsplitDMD(data,z(:,p(1,k):p(1,k+1)),s_1,t_1,iter,r,t);
                t_split = [ t_split t_1 ];
                s_k = s_k + s_1 - 1;
                iter = iter + iter1;
            else
                return;
            end
            if k == (s-1) && abs(p(1,k+2) - p(1,k+1)) > del
                stall = false;
                t_2 = [ t_split round((p(1,k+2) + p(1,k+1))/2) ];
                [s_2,t_2,iter2] = run_nsplitDMD(data,z(:,p(1,k+1):p(1,k+2)),s_2,t_2,iter,r,t);
                t_split = [ t_split t_2 ];
                s_k = s_k + s_2 - 1;
                iter = iter + iter2;
            else
                return;
            end
            % otherwise, interval is too short to split
        end
    end

    if stall
        iter = iter + 1;
        s = s + 1;
        t_split = zeros(1,s);
        [s_k,t_split,iter3] = run_nsplitDMD(data,z,s,t_split,iter,r,t);
        iter = iter + iter3;
    end

    t_split = sort(t_split);
    t_split = unique(t_split.','rows').';
    n_split = size(t_split,2);  
    
    test_t = size(t_split,2)-2;
    
    display(t_split)
    if size(t_split,2) > 2
        ct = 0;
        t_split1 = t_split;
        for k = 1:test_t
            if t_split(1,k+1) == (t_split(1,k) + 1)
                t_split1 = [t_split1(1,1:k-ct) t_split1(1,(k+2-ct):end)];
                ct = ct + 1;
                n_split = n_split - 1;
            end
        end
        t_split = t_split1;
    end
        
    if (t_split(1,end-1) + 1) == t_split(1,end)
        t_split = t_split(1,1:end-1);
        n_split = n_split - 1;
    end

      
    
end