function [STC, AUC, MCC, F1S, OLP, estTh]  = anlys_STC(curryloc, currycdr, source, Phi, stc_option, disp_opt)
% thin function statistically analyzes the solution
Th  = linspace(1, 0, 1000);

%%
[~, simind]  = find_nvoxel(ppatches(source.epl, 'row'), curryloc); 
simind       = ismember(1:size(curryloc, 2), unique(simind)).';


%%
switch stc_option
    
    case 'Norm'
        cdr  = norms(squeeze(currycdr(4, :, :)), 2, 2);
    
    case 'MGFP'
        [P, TP]  = findpeaks(std(Phi, 1, 1)); 
        TP       = TP(P >= max(P)/sqrt(10));
        for t = 1 : length(TP)
            cdr(:, t)  = abs(squeeze(currycdr(4, :, TP(t)))).';
        end

end


%%
for t = 1 : size(cdr, 2)
    
    for i_th = 1 : length(Th)

        [TP(i_th, t), TN(i_th, t), FP(i_th, t), FN(i_th, t)] = find_rates(cdr(:, t), simind, Th(i_th));

    end
    
end

[AUC, MCC, F1S, OLP, estTh] = find_STC(Th, TP, TN, FP, FN, disp_opt);

STC.Th  = Th;
STC.TP  = TP;
STC.TN  = TN;
STC.FP  = FP;
STC.FN  = FN;

    
end


function [TP, TN, FP, FN] = find_rates(map, simind, Thr)
    % refer to https://en.wikipedia.org/wiki/Precision_and_recall and
    % https://en.wikipedia.org/wiki/Matthews_correlation_coefficient

    %=====================================================================%
    estind  = (map >= Thr*max(map));
    TP  = sum((simind == 1) & (estind==1));
    TN  = sum((simind == 0) & (estind==0));
    FP  = sum((simind == 0) & (estind==1));
    FN  = sum((simind == 1) & (estind==0));
    
end


