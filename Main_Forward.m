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

% This script generates and saves electro-potential maps by solving the 
% forward problem for a set of source configurations, leadfield matrices,
% additive noise and for (different/a single) case(s) in nested for-loops. 
% 
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

clc, close all, clear all
%%

cfg  = [1, 2, 3, 4]; % config
lfd  = [29]; % leadfield
noi  = [3, 9]; % noise
src  = [1]; % source

n_cfg  = length(cfg);
n_lfd  = length(lfd);
n_noi  = length(noi);
n_src  = length(src);

%=========================================================================%%=========================================================================%
h = waitbar(0, 'Preparing to solve forward problem...');
for i_cfg = 1 : n_cfg

    Main_Setup
    start_point  = tic;
    BaselineLen  = 1000;
    load('Resource/SOURCE_Example.mat')
    
    for i_lfd = 1 : n_lfd

        load(['Resource/' lfd_forward{lfd(i_lfd), 2} '.mat'])

        for i_noi = 1 : n_noi
            
            rng(1)

            for i_src = 1 : n_src
                
                %=========================================================%%=========================================================%
                Phi  = solve_forward_problem(currylfd, curryloc, SOURCE(cfg(i_cfg), lfd(i_lfd), src(i_src)), lfd_forward{lfd(i_lfd), 4});
                
                %=========================================================%%=========================================================%
                [Phi, ~]  = add_noise(Phi, noise_setup{noi(i_noi), 2}, noise_setup{noi(i_noi), 3}, BaselineLen, currylfd, [], 0.5, ...
                    SOURCE(cfg(i_cfg), lfd(i_lfd), src(i_src)).acr);
                save(['Resource/DAT Files/phi_' num2str(cfg(i_cfg)) '_' lfd_forward{lfd(i_lfd), 3} '_' ...
                    noise_setup{noi(i_noi), 4} '_' num2str(src(i_src)) '_Example.dat'], 'Phi', '-tabs', '-ascii')
                
                %=========================================================%%=========================================================%
                percentage  = ((i_lfd-1)*n_noi*n_src + (i_noi-1)*n_src + i_src) / (n_lfd*n_noi*n_src);
                waitbar( percentage, h, ['Solving forward problem (' num2str(cfg(i_cfg)) '), ' datestr((1-percentage)/percentage*toc(start_point)/60/60/24, 'HH:MM:SS') ' left...']); 

            end

        end

    end
    
end
