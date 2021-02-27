function [AUC, OLP, Th_OLP]  = anlys_AUC(curryloc, currycdr, source, Phi, auc_option, disp_opt)

% Th  = linspace(1, 0, 1000);
Th  = logspace(log10(1), log10(1e-6), 1000);
[~, simind] = find_nvoxel(ppatches(source.epl, 'row'), curryloc); 
simind  = ismember(1:size(curryloc, 2), unique(simind)).';



switch auc_option
    
    case 'MGFP'
        [P, T]  = findpeaks(std(Phi, 1, 1)); 
        T       = T(P >= max(P)/sqrt(10));
        TPR     = zeros(length(Th), length(T));
        FPR     = zeros(length(Th), length(T));
        AUC     = zeros(length(T), 1);
        OLP     = zeros(length(T), 3);
        Th_OLP  = zeros(length(T), 1);

        for t = 1 : length(T)
            cdr  = abs(squeeze(currycdr(4, :, T(t)))).';
            for i_Th = 1 : length(Th)

                [TPR(i_Th, t), FPR(i_Th, t), ~] = find_rates(cdr, simind, Th(i_Th));

            end
            AUC(t)      = trapz(FPR(:, t), TPR(:, t));
            [Th_OLP(t), FPR_C(t), ~]   = l_corner(FPR(:, t), 1-TPR(:, t), Th);
            [~, ~, olp] = find_rates(cdr, simind, Th_OLP(t));
            OLP(t, :)   = [olp 2*olp(1)*olp(2)/(olp(1)+olp(2))];

        end
        
    case 'Norm'
        TPR = zeros(length(Th), 1);
        FPR = zeros(length(Th), 1);
        cdr  = norms(squeeze(currycdr(4, :, :)), 2, 2);
        for i_Th = 1 : length(Th)
            [TPR(i_Th), FPR(i_Th), ~] = find_rates(cdr, simind, Th(i_Th));
        end
        AUC    = trapz(FPR, TPR);
        [Th_OLP, FPR_C, ~]  = l_corner(FPR, 1-TPR, Th);
        [~, ~, olp] = find_rates(cdr, simind, Th_OLP);
        OLP  = [olp 2*olp(1)*olp(2)/(olp(1)+olp(2))];

end


if(strcmp(disp_opt, 'Display'))
    figure, 
    plot(FPR, TPR)
    T  = length(Th_OLP);
    hold on, line([FPR_C.'; FPR_C.'], [zeros(1, T); ones(1, T)])
    ylabel('True Positive Rate')
    xlabel('False Positive Rate')
    for i = 1 : T
        lgnd{i} = ['Th = ' num2str(Th_OLP(i))];
    end
    legend(lgnd)
end
Th_OLP
    
end


function [TPR, FPR, OLP] = find_rates(map, simind, Th)
    % refer to https://en.wikipedia.org/wiki/Precision_and_recall

    estind  = (map >= Th*max(map));
    
    TP  = sum((simind == 1) & (estind==1));
    TN  = sum((simind == 0) & (estind==0));
    FP  = sum((simind == 0) & (estind==1));
    FN  = sum((simind == 1) & (estind==0));
    
    TPR  = TP / (TP + FN); % sensitivity
    FPR  = FP / (FP + TN); % specificity complemntary
    
    OLP(1)  = TP / (TP + FN); % area overlap over true area or "recall"
    OLP(2)  = TP / (TP + FP); % area overlap over estimated area or "precision:
    
end