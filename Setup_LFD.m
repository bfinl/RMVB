% This script sets up the uncertainty leadfield matrices to estimate the 
% uncertainty regions (ellipsoids) as well as the leadfield matrix for 
% solving the inverse problem.
% 
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

%%
clc, close all, clear all,

%% #1    

    rng(1)
    forwardlfd = [];
    load('Resource/curry_lfd_BEM-0.0400_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0417_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0435_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0455_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0476_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0500_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0526_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0556_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0588_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0625_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0667_CTX-1.1mm-fxd_EEG-Biosemi-128.mat');
    forwardlfd  = cat(3, forwardlfd, revise_lfd_dir(currylfd, 'fxd'));  
    forwardloc  = currylfd(1:3, :);
    forwardori  = currylfd(4:6, :);

    inverselfd = [];
    load('Resource/curry_lfd_BEM-0.0400_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0417_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0435_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0455_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0476_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0500_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0526_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0556_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0588_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0625_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    load('Resource/curry_lfd_BEM-0.0667_CTX-5.0mm-fxd_EEG-Biosemi-128.mat');
    inverselfd  = cat(3, inverselfd, revise_lfd_dir(currylfd, 'fxd'));
    inverseloc  = currylfd(1:3, :);
    inverseori  = currylfd(4:6, :);
    
    m  = size(currylfd, 2); 
    if(norm(curryloc(:, 1:m)-currylfd(1:3, :), 'fro'))
        display('Warning: the off inner skull grids are not expanded at the end.')
    end
    dsplyloc  = curryloc(:, m+1:end);
    curryori  = currylfd(4:6, :);
    curryloc  = currylfd(1:3, :);
    curryani  = [];
    load('Resource/SRCIDXMAT-fxd.mat')
    
    %=====================================================================%
    % This part is to deal with the coarse griding of the invese problem
    % which makes the orientations inconsistent over the neibourhood.
    [~, IDX]  = find_nvoxel(inverseloc, forwardloc);
    curryori  = forwardori(:, IDX);
    curryloc  = forwardloc(:, IDX);
    inverseori  = forwardori(:, IDX);
    inverseloc  = forwardloc(:, IDX);
    inverselfd  = forwardlfd(:, IDX, :);
    %=====================================================================%
        
%     [gammamat, currylfd]  = find_uncertainty(forwardlfd, inverselfd, forwardloc, inverseloc, forwardori, inverseori, 'TEST', 'FLD', 'fxd', 1:m);
%     orthglfd  = find_orthg_lfd(currylfd, 'fxd');
% % % %     [SDI, TAG]  = generate_seed_points(currylfd, currylfd, orthglfd, currytri, 4, 512, 'GROVA', 'fxd');
% % % %     [curryani, PMP]  = generate_parcellation(SDI, TAG, currytri, size(curryloc, 2));
% % % %     for i = 1 : length(SDI),  SDL{i}  = curryloc(:, SDI(i)); end
% % % %     setup_display(PMP, currytri, curryloc, dsplyloc, 0, '', SDL, 'pnt');
%     save('Resource/lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128_TEST.mat', 'currylfd', 'curryloc', 'curryori', 'currytri', 'dsplyloc', 'orthglfd', 'gammamat', 'curryani');
    
%     [gammamat, currylfd]  = find_uncertainty(forwardlfd, inverselfd, forwardloc, inverseloc, forwardori, inverseori, 'EMP', 'FLD', 'fxd', 1:m);
%     orthglfd  = find_orthg_lfd(currylfd, 'fxd');
% %     [SDI, TAG]  = generate_seed_points(currylfd, currylfd, orthglfd, currytri, 4, 256, 'GROVA', 'fxd');
% %     [curryani, PMP]  = generate_parcellation(SDI, TAG, currytri, size(curryloc, 2));
% %     for i = 1 : length(SDI), SDL{i}  = curryloc(:, SDI(i)); end
% %     setup_display(PMP, currytri, curryloc, dsplyloc, 0, '', SDL, 'pnt');
%     save('Resource/lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128_EMP.mat', 'currylfd', 'curryloc', 'curryori', 'currytri', 'dsplyloc', 'orthglfd', 'gammamat', 'curryani');
    
    [gammamat, currylfd]  = find_uncertainty(forwardlfd, inverselfd, forwardloc, inverseloc, forwardori, inverseori, 'CAM', 'FLD', 'fxd', 1:m);
    orthglfd  = find_orthg_lfd(currylfd, 'fxd');
% %     [SDI, TAG]  = generate_seed_points(currylfd, currylfd, orthglfd, currytri, 4, 512, 'GROVA', 'fxd');
% %     [curryani, PMP]  = generate_parcellation(SDI, TAG, currytri, size(curryloc, 2));
% %     for i = 1 : length(SDI),  SDL{i}  = curryloc(:, SDI(i)); end
% %     setup_display(PMP, currytri, curryloc, dsplyloc, 0, '', SDL, 'pnt');
    save('Resource/lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128_CAM.mat', 'currylfd', 'curryloc', 'curryori', 'currytri', 'dsplyloc', 'orthglfd', 'gammamat', 'curryani');

%     [gammamat, currylfd]  = find_uncertainty(forwardlfd, inverselfd, forwardloc, inverseloc, forwardori, inverseori, 'AMIR', 'FLD', 'fxd', 1:m);
%     orthglfd  = find_orthg_lfd(currylfd, 'fxd');
% %     [SDI, TAG]  = generate_seed_points(currylfd, currylfd, orthglfd, currytri, 4, 512, 'GROVA', 'fxd');
% %     [curryani, PMP]  = generate_parcellation(SDI, TAG, currytri, size(curryloc, 2));
% %     for i = 1 : length(SDI),  SDL{i}  = curryloc(:, SDI(i)); end
% %     setup_display(PMP, currytri, curryloc, dsplyloc, 0, '', SDL, 'pnt');
%     save('Resource/lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128_AMIR.mat', 'currylfd', 'curryloc', 'curryori', 'currytri', 'dsplyloc', 'orthglfd', 'gammamat', 'curryani');


  