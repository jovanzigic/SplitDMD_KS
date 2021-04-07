%100 points is 500 steps
tn = t;
for k = 1:50000
    tn(end+1) = tn(end) + 0.02;
end

%1. integrate Phi over space
%2. 

% w1 = cell2mat(w(1,10));
% e1 = cell2mat(e(1,10));
% b1 = cell2mat(b(1,10));

% Xd1n_1 = w1*diag(b1)*exp(e1*tn(1:12499));
% Xd1n_3 = w1*diag(b1)*exp(e1*tn(20001:25001));
% Xd1n_2 = w2*diag(b2)*exp(e2*tn(12500:20000));
% % Xd1n_4 = w2*diag(b2)*exp(e2*tn(32500:40001));
% Xd1n = [Xd1n_1 Xd1n_2 Xd1n_3 ];

Xd1n_1 = w*diag(b)*exp(e*tn);
% Xd1n_3 = w1*diag(b1)*exp(e1*tn(20001:25001));
% Xd1n_2 = w2*diag(b2)*exp(e2*tn(12500:20000));
% Xd1n_4 = w2*diag(b2)*exp(e2*tn(32500:40001));
Xd1n = [Xd1n_1 ];

realXd1n = real(Xd1n);
% realXd1n = realXod1 - realXod;

Nx = size(realXd1n,1);
Nt = size(realXd1n,2);

xp = linspace(0,1,Nx);
% tp = linspace(0,interval,Nt);

figure
% mesh(x,tn(19001:25001),realXd1n'); 
mesh(x,tn,realXd1n'); 
%view([.7 -.8 .6])
view([0 0 1])
title(['Optimized DMD Model with Prediction, order = ' , num2str(r_dim) ] )
xlabel('x'); ylabel('t'); zlabel('z_{OD}')
xlim([0 1]) 
% ylim([0 420])