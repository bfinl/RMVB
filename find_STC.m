function [auc, mcc, f1s, olp, estTh] = find_STC(Th, TP, TN, FP, FN, disp_opt)


%=========================================================================%
PPV  = TP ./ (TP  + FP); % precision/overlap normalized by estimated area
TPR  = TP ./ (TP  + FN); % sensitivity/overlap normalized by true area
FPR  = FP ./ (FP  + TN); % specificity complemntary

% Matthews correlation coefficient
%=========================================================================%
MCC  = (TP.*TN - FP.*FN)./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
MCC(isinf(MCC) | isnan(MCC))  = 0; 

% F1 scores
%=========================================================================%
F1S  = 2*PPV.*TPR./(PPV+TPR);

      
%%
[mcc, imcc]  = max(MCC, [], 1);
[f1s, if1s]  = max(F1S, [], 1); 
estTh        = Th(imcc);

for t = 1 : size(TPR, 2)
    auc(t)     = trapz(FPR(:, t), TPR(:, t));
    olp(t, :)  = [TPR(imcc(t), t) PPV(imcc(t), t)];
    dFPR(t)    = FPR(imcc(t), t);
end


%%
if(strcmp(disp_opt, 'Display'))
    figure, 
    plot(FPR, TPR)
    TP  = length(estTh);
    hold on, line([dFPR; dFPR], [zeros(1, TP); ones(1, TP)])
    ylabel('True Positive Rate')
    xlabel('False Positive Rate')
    for i = 1 : TP
        lgnd{i} = ['Th = ' num2str(estTh(i))];
    end
    legend(lgnd)
end

end