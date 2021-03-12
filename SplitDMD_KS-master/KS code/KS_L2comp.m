clear, close all

interval = 400; % adjust parameter
re = 402.3; % adjust parameter
r_dim = 13; % adjust parameter

split = 1;

if split == 1
    fileload = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split2_2');
else
    fileload = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
end
load(fileload)

% Simulate DMD ROM and compute L2 error
[errorL2_DMD,errorL2_OD] = sim_DMD_KS(realXdmd,errorDMD,realXod,errorOD,r_dim,interval,M);

% Save data
if split == 1
    filesave = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim),'_split2_2');
else
    filesave = strcat('ROM_KS_data_t',int2str(interval),'_L',int2str(100*re),'_r',int2str(r_dim));
end
save(filesave)