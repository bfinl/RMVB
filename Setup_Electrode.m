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

% This script sets up electrodes based on a Biosemi cap cofiguration with 128 channels.
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

%%
clc, close all, clear all,

%% 
ELECTRODE = struct('loc', [], 'spc', '');


%% Electrode Configuration # 1
load('Resource\Biosemi-128.mat')
electrode.spc = 'Biosemi-128';
electrode.loc = curryloc;
ELECTRODE(1)  = electrode;
elec_save   = electrode.loc.';
save(['Resource\' electrode.spc '.txt'], 'elec_save', '-ASCII');
figure, scatter3(electrode.loc(1, :), electrode.loc(2, :), electrode.loc(3, :), 'o')


%%
save('Resource\ELECTRODE', 'ELECTRODE')