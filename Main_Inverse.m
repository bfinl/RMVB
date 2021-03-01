% This program is a free software for academic research: you can redistribute 
% it and/or modify it for non-commercial uses, under the license terms provided 
% with the package at the GitHub page where this package is downloaded. 
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the License for more details. This program is 
% for research purposes only. This program CAN NOT be used for commercial purposes. 
% This program SHOULD NOT be used for medical purposes. The authors WILL NOT 
% be responsible for using the program in medical conditions.
% ==========================================

% This script reconstructs and saves source activities by solving the inverse 
% problem using different reconstruction methods for a set of source configurations, 
% leadfield matrices, additive noise and for (different/a single) 
% case(s) in nested for-loops. 
% 
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

% clc, close all, clear all
%%

cfg  = [4]; % 4: extneded activity with 3 nodes (see setup_source_vExample.m)
lfd  = [29]; % 29: the actual leadfield (used in naming files)(see main_setup.m)
noi  = [9]; % 9: 20dB AWGN (see main_setup.m)
src  = [1]; % 1: first source activity generated. This is 1 for the example version. (see setup_source_vExample.m)
cdr  = {'LCMV', 'LCMV-ND-DN', 'RMVB', 'RMVB-ND-DN'}; % reconstruction methods
dle_opt  = {'PCA-FOC', 'PCA-FOC', 'PCA-FOC', 'PCA-EXT'}; % Focal or extended PCA application. For cfg=4, 'PCA-EXT' is used. 
ThL  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; % thresholing value

cdr_disp_opt  = 1;
FilterOption  = 'FiltFilt';

n_cfg  = length(cfg);
n_lfd  = length(lfd);
n_noi  = length(noi);
n_src  = length(src);
n_cdr  = length(cdr);

elapsed_time  = zeros(n_cfg, 1);
RLTS  = cell(n_cfg, n_lfd, n_noi, n_src, n_cdr, 6);
DUR  = NaN;
CNN  = NaN;


%%  Main
h = waitbar(0, 'Preparing to solve inverse problem...');
for i_cfg = 1 : n_cfg
    
    rng(1)
    Main_Setup
    cvx_solver sedumi
    BaselineLen = 1000;
    start_point  = tic;
    addpath('glasso-library')
    addpath('regu-customized')
    addpath('NewtonRaphson_0.5')
    load('Resource/ELECTRODE.mat')
    load('Resource/SOURCE_Example.mat')

    for i_lfd = 1 : n_lfd
        
%         if(lfd(i_lfd)==26 && cfg(i_cfg)~=4), continue, end
        if(lfd(i_lfd)==36 && cfg(i_cfg)~=2), continue, end
        load(['Resource/' lfd_inverse{lfd(i_lfd), 2} '.mat']); 

        for i_noi = 1 : n_noi

            for i_src = 1 : n_src

                source  = SOURCE(cfg(i_cfg), lfd(i_lfd), src(i_src));
                Phi  = load(['Resource/DAT Files/phi_' num2str(cfg(i_cfg)) '_' lfd_forward{lfd(i_lfd), 3} '_' noise_setup{noi(i_noi), 4} '_' num2str(src(i_src)) '_Example.dat']);
                Phi  = my_filter(Phi, source.frq, [1, 50], FilterOption, BaselineLen, 'No Display');
                Phi_az  = Phi(:, SOURCE(cfg(i_cfg), lfd(i_lfd), src(i_src)).acr+BaselineLen);
                Phi_bl  = Phi(:, 1:BaselineLen);

                for i_cdr = 1 : n_cdr
                                 
                    %===================================================================================================================================================================%
                    tic, [CURRYCDR] = solve_inverse_problem(Phi_bl, Phi_az, currylfd, orthglfd, curryori, currytri, curryani, lfd_inverse{lfd(i_lfd), 4}, cdr{i_cdr}, source.nnd, ...
                        source, curryloc, dsplyloc, gammamat); DUR  = toc;
                     
                    lbl  = [cdr{i_cdr} '_' num2str(cfg(i_cfg)) '_' lfd_forward{lfd(i_lfd), 3} '_' lfd_inverse{lfd(i_lfd), 3} '_' noise_setup{noi(i_noi), 4} '_' num2str(src(i_src))];
                    currycdr = CURRYCDR;
                    save(['Resource/DAT Files/cdr_' lbl '_Example.mat'], 'currycdr', 'curryloc')
                    load(['Resource/DAT Files/cdr_' lbl '_Example.mat'])

                    %===================================================================================================================================================================%
                    [STC, AUC, MCC, F1S, OLP, estTh]  = anlys_STC(curryloc, currycdr, source, Phi_az, 'Norm', 'No Display');

                    %===================================================================================================================================================================%
                    [SDL, EPL, DLE, DLI]  = anlys_DLE(curryloc, currycdr, currytri, source, ThL(i_cdr), mean(estTh), 'Predefined', dle_opt{cfg(i_cfg)});

                    %===================================================================================================================================================================%
                    [CRR1, SNR1]  = anlys_TCR(source.epl, source, currycdr, curryloc, 1:source.nnd, 'svd-v2', 'Normalize', 'No Display');
                    [CRR2, SNR2]  = anlys_TCR(EPL, source, currycdr, curryloc, DLI, 'svd-v2', 'Normalize', 'No Display');
                    TCR.CRR1  = CRR1; TCR.CRR2  = CRR2; TCR.SNR1  = SNR1; TCR.SNR2  = SNR2;

                    %===================================================================================================================================================================%
%                     [COV]  = anlys_COV(source.epl, EPL, source.epc, currycdr, curryloc);

                    %===================================================================================================================================================================%
%                         [KLD]  = anlys_KLD(currycdr, curryloc, source, Phi_az);

                    %===================================================================================================================================================================%
%                     [mnt] = find_accurate_cdr(SDL, EPL, Phi_bl, Phi_az, curryloc, currylfd, lfd_inverse{lfd(i_lfd), 4}, 'ext', 'Display', source.frq, DLI, currycdr);
%                     [gm2, CNN]  = anlys_CNN(mnt, 10, 40, 1, 50, source.frq, 1000, 0.05, 'Normalize', 'fpe', 'Display', source.dtf(DLI, DLI, :));

                    %===================================================================================================================================================================%
                    if(cdr_disp_opt) 
%                             setup_display(currycdr, currytri, curryloc, dsplyloc, mean(estTh), lbl, SDL, 'pnt');
%                         setup_display(currycdr, currytri, curryloc, dsplyloc, mean(estTh), 'Recovered Patches!', EPL, 'net');
                        setup_display(currycdr, currytri, curryloc, dsplyloc, mean(estTh), lbl, source.epl, 'net');
                        setup_display(currycdr, currytri, curryloc, dsplyloc, 0.01, lbl, source.epl, 'net');
                    end

                    %===================================================================================================================================================================%
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 1}  = DUR;
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 2}  = DLI;
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 3}  = DLE;  
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 4}  = STC;
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 5}  = TCR;
                    RLTS{i_cfg, i_lfd, i_noi, i_src, i_cdr , 6}  = CNN;
                    fprintf('CDR: %s\nDUR: %s\nDLE: %s\nAUC: %s\nMCC: %s\nF1S: %s\nSNR1: %s\nSNR2: %s\n',...
                        cdr{i_cdr}, my_vec2str(DUR), my_vec2str(DLE), my_vec2str(AUC), my_vec2str(MCC), my_vec2str(F1S), my_vec2str(SNR1), my_vec2str(SNR2))
                    display('____________________________________________________________________________________________')

                    %===================================================================================================================================================================%
                    percentage = ((i_lfd-1)*n_noi*n_src*n_cdr + (i_noi-1)*n_src*n_cdr + (i_src-1)*n_cdr + i_cdr) / (n_lfd*n_noi*n_src*n_cdr);
                    waitbar(percentage, h, ['Inverse (' num2str(cfg(i_cfg)) ')... ' datestr(toc(start_point)/60/60/24, 'HH:MM:SS') ' done... ' datestr((1-percentage)/percentage*toc(start_point)/60/60/24, 'HH:MM:SS') ' left...']); 

                end

            end

        end       

    end
    
    elapsed_time(cfg(i_cfg))  = toc(start_point);
    
end

close(h)
fprintf('Elapsed Time: %s\n', my_vec2str(elapsed_time))
