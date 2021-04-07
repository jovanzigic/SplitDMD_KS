clear, close all

interval = 400; % adjust parameter
re = 12.6; % adjust parameter
r_dim = 11; % adjust parameter

split = 1;

if split == 1
    fileload1 = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split');
else
    fileload1 = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
end
load(fileload1)

% Simulate DMD ROM and compute L2 error
[errorL2_DMD,errorL2_OD] = sim_DMD_KS(realXdmd,errorDMD,realXod,errorOD,r_dim,interval,M);

% Save data
save(fileload1)