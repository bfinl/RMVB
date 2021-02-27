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