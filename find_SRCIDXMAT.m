Main_Setup
load('Resource\SOURCE.mat')
load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128_CAM.mat');
SRCIDXMAT = [];
cfg  = [1, 5, 9];
lfd  = [23 : 35];
temp1 = [];
for i_cfg = 1 : length(cfg)
    for i_lfd = 1 : length(lfd)
        for i_src = 1 : 50
            [~, ind]  = find_nvoxel(ppatches(SOURCE(cfg(i_cfg), lfd(i_lfd), i_src).epl, 'row'), curryloc);
            SRCIDXMAT = [SRCIDXMAT ind];
            temp1 = [temp1 ind];
        end
    end
end
SRCIDXMAT = unique(SRCIDXMAT);
temp1 = unique(temp1);

save('Resource\SRCIDXMAT-fxd.mat', 'SRCIDXMAT')

%%
Main_Setup
load('Resource\SOURCE.mat')
load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-rot_EEG-Biosemi-128_CAM.mat');
SRCIDXMAT = [];
cfg  = [1, 5, 9];
lfd  = [1 : 22];
temp1 = [];
for i_cfg = 1 : length(cfg)
    for i_lfd = 1 : length(lfd)
        for i_src = 1 : 50
            [~, ind]  = find_nvoxel(ppatches(SOURCE(cfg(i_cfg), lfd(i_lfd), i_src).epl, 'row'), curryloc);
            SRCIDXMAT = [SRCIDXMAT ind];
            temp1 = [temp1 ind];
        end
    end
end
SRCIDXMAT = unique(SRCIDXMAT);
temp1 = unique(temp1);

save('Resource\SRCIDXMAT-rot.mat', 'SRCIDXMAT')