% This script sets up 4 source activity configurations in the brain 
% according to section "Computer Simulation Protocol" of the RMVM paper
% (DOI: 10.1109/TBME.2018.2859204)
% 
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

%%
clc, close all, clear all,

%%
Main_Setup
h = waitbar(0, 'Waiting');
lfd  = 1 : size(lfd_forward, 1);
lfd  = [29];

source.frq  = [];   % frequency
source.nnd  = [];   % number of nodes in the network
source.nts  = [];   % number of time series points
source.sdl  = [];   % seed location
source.sdo  = [];   % seed orientation
source.sdm  = [];   % seed moment
source.epl  = [];   % entire patch location
source.epo  = [];   % entire patch orientation
source.epm  = [];   % entire patch moment
source.epa  = [];   % entire patch area
source.dtf  = [];   % dtf connectivity
source.acr  = [];   % active region
source.spc  = '';   % specificaion
            
                                                                                                                      
for i_lfd = 1 : length(lfd)

    
    start_point  = tic;
    load(['Resource/' lfd_forward{lfd(i_lfd), 2} '.mat'])
    if(i_lfd == 1), disp_opt = 'No Display'; else disp_opt = 'No Display'; end;
    
    
    %% Config. # 1
    rng(1)
    cfg  = 1;
    frq  = 1000;
    nts  = 0.4 * frq;
    nnd  = 1;
    sdm  = setup_activity(frq, nnd, nts, 'Regular1-v2', 'Normalize', disp_opt);
    dtf  = [];
    acr  = 1:nts;
    spc  = 'Regular1-v2';
    stp  = tic;
    
    rng(1)
    loc  = [74.3,   15.4,  38.7].';
    [~, ind]  = find_nvoxel(loc, curryloc); 
    [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch(curryloc, curryori, currytri, nnd, sdm, 'area', [0], ind);

    source.frq  = frq; 
    source.nnd  = nnd; 
    source.nts  = nts;   
    source.sdl  = sdl; 
    source.sdo  = sdo;
    source.sdm  = sdm;
    source.epl  = epl; 
    source.epo  = epo;
    source.epm  = [];
    source.epa  = epa;
    source.dtf  = dtf;
    source.epc  = epc;
    source.acr  = acr;
    source.spc  = spc;
    
    SOURCE(cfg, lfd(i_lfd), 1)  = source;
 
    
    %% Config. # 2
    rng(1)
    cfg  = 2;
    frq  = 1000;  
    nts  = 0.4 * frq;
    nnd  = 3;
    sdm  = setup_activity(frq, nnd, nts, 'Regular3-v4', 'Normalize', disp_opt);
    dtf  = [];
    acr  = 1:nts;
    spc  = 'Regular3-v4';
    stp  = tic;

    rng(1)
    loc  = [-10.5,  80.2,  99.8;  
             74.3,  15.4,  38.7;
            -14.1, -55.6,  72.7;
            ].';
    [~, ind]  = find_nvoxel(loc, curryloc); 
    [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch(curryloc, curryori, currytri, nnd, sdm, 'area', [0, 0, 0], ind);

    source.frq  = frq; 
    source.nnd  = nnd; 
    source.nts  = nts;   
    source.sdl  = sdl; 
    source.sdo  = sdo;
    source.sdm  = sdm;
    source.epl  = epl; 
    source.epo  = epo;
    source.epm  = [];
    source.epa  = epa;
    source.dtf  = dtf;
    source.epc  = epc;
    source.acr  = acr;
    source.spc  = spc;
    
    SOURCE(cfg, lfd(i_lfd), 1)  = source;
    
    
    %% Config. # 3
    rng(1)
    cfg  = 3;
    frq  = 1000;  
    nts  = 0.4 * frq;
    nnd  = 3;
    sdm  = setup_activity(frq, nnd, nts, 'Regular3-v5', 'Normalize', disp_opt);
    dtf  = [];
    acr  = 1:nts;
    spc  = 'Regular3-v5';
    stp  = tic;

    rng(1)
    loc  = [-10.5,  80.2,  99.8;  
             74.3,  15.4,  38.7;
            -14.1, -55.6,  72.7;
            ].';
    [~, ind]  = find_nvoxel(loc, curryloc); 
    [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch(curryloc, curryori, currytri, nnd, sdm, 'area', [0, 0, 0], ind);

    source.frq  = frq; 
    source.nnd  = nnd; 
    source.nts  = nts;   
    source.sdl  = sdl; 
    source.sdo  = sdo;
    source.sdm  = sdm;
    source.epl  = epl; 
    source.epo  = epo;
    source.epm  = [];
    source.epa  = epa;
    source.dtf  = dtf;
    source.epc  = epc;
    source.acr  = acr;
    source.spc  = spc;
    
    SOURCE(cfg, lfd(i_lfd), 1)  = source;


    %% Config. # 4
    rng(1)
    cfg  = 4;
    frq  = 1000;  
    nts  = 0.4 * frq;
    nnd  = 3;
    sdm  = setup_activity(frq, nnd, nts, 'Regular3-v4', 'Normalize', disp_opt);
    dtf  = [];
    acr  = 1:nts;
    spc  = 'Regular3-v4';
    stp  = tic;

    rng(1)
    loc  = [-10.5,  80.2,  99.8; 
             74.3,  15.4,  38.7;
            -14.1, -55.6,  72.7;
            ].';
    [~, ind]  = find_nvoxel(loc, curryloc); 
    [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch(curryloc, curryori, currytri, nnd, sdm, 'area', [750, 750, 750], ind);

    source.frq  = frq; 
    source.nnd  = nnd; 
    source.nts  = nts;   
    source.sdl  = sdl; 
    source.sdo  = sdo;
    source.sdm  = sdm;
    source.epl  = epl; 
    source.epo  = epo;
    source.epm  = [];
    source.epa  = epa;
    source.dtf  = dtf;
    source.epc  = epc;
    source.acr  = acr;
    source.spc  = spc;
        
    SOURCE(cfg, lfd(i_lfd), 1)  = source;
    
    
%%

    elapsed_time(lfd(i_lfd))  = toc(start_point);
    
    
end


%%
save('Resource/SOURCE_Example', 'SOURCE', '-v7.3')
