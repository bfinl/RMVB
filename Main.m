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

clc, close all, clear all,
%%
Setup_Electrode

%%
Setup_LFD

%%
Setup_Source

%%
Main_Forward
Main_Inverse
save('Results', '-v7.3')








