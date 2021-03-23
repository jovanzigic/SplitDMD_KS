function [n_split,t_split,iter] = run_nsplit(data,z,s,t_split,iter)
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

maxiter = 20;

n_split = s;
    
if iter < maxiter

    % choose random space across time
    x_rand = randi(size(z,1),1,1);

    x = x_rand;

    eps = round((1/5)*range(data(x,:)));
    del = round((1/5)*size(data(x,:),2));

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

    s_1 = 1; 
    s_2 = 1; 
    s_k = s;    

    for k = 1:(s-1) % usually, s = 2
        if (abs( max(z(x,p(1,k+1):p(1,k+2))) - max(z(x,p(1,k):p(1,k+1))) ) > eps ...
                || abs( min(z(x,p(1,k+1):p(1,k+2))) - min(z(x,p(1,k):p(1,k+1))) ) > eps)
            if abs(p(1,k+1) - p(1,k)) > del
                stall = false;
                t_1 = [ t_split round((p(1,k+1) + p(1,k))/2) ];
                [s_1,t_1,iter1] = run_nsplit(data,z(:,p(1,k):p(1,k+1)),s_1,t_1,iter);
                t_split = [ t_split t_1 ];
                s_k = s_k + s_1 - 1;
                iter = iter + iter1;
            else
                return;
            end
            if k == (s-1) && abs(p(1,k+2) - p(1,k+1)) > del
                stall = false;
                t_2 = [ t_split round((p(1,k+2) + p(1,k+1))/2) ];
                [s_2,t_2,iter2] = run_nsplit(data,z(:,p(1,k+1):p(1,k+2)),s_2,t_2,iter);
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
        s_j = s + 1;
        t_split = zeros(1,s_j);
        [~,t_split,iter3] = run_nsplit(data,z,s_j,t_split,iter);
        iter = iter + iter3;
    end

    t_split = sort(t_split);
    t_split = unique(t_split.','rows').';
    n_split = size(t_split,2);  
    
    test_t = size(t_split,2)-2;
    
%     display(t_split)
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
        
        if (t_split(1,end-1) + 1) == t_split(1,end)
            t_split = t_split(1,1:end-1);
            n_split = n_split - 1;
        end
    end
    
end

end