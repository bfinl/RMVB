function [CRR, SNR]  = anlys_TCR(epl, source, currycdr, curryloc, dli, ori_revision_opt, nlz_opt, dsply_opt)

frq  = source.frq;
sdm  = source.sdm(dli, :);
mnt  = zeros(size(sdm));
nnd  = size(sdm, 1);
nts  = size(sdm, 2);
mnt_cmp  = zeros([size(sdm), 3]);

for i = 1 : length(epl)
    [~, ind]  = find_nvoxel(epl{i}, curryloc);
    tmprycdr  = revise_cdr_ori(currycdr(:, ind, :), ori_revision_opt);
    mnt(i, :) = squeeze(mean(tmprycdr(4, :, :), 2));
    
    for j = 1 : 3
        mnt_cmp(i, :, j) = squeeze(mean(currycdr(4, ind, :).*currycdr(j, ind, :), 2));
    end
end

%%
% for i = 1 : length(epl)
%     [~, ind]  = find_nvoxel(epl{i}, curryloc);
%     tmprycdr  = currycdr(:, ind, :);
%     mnt(i, :) = squeeze(abs(mean(tmprycdr(4, :, :), 2)));
%     sdm(i, :) = abs(sdm(i, :));
%     
%     for j = 1 : 3
%         mnt_cmp(i, :, j) = squeeze(mean(currycdr(4, ind, :).*currycdr(j, ind, :), 2));
%     end
% end


%%

if(strcmp(nlz_opt, 'Normalize'))
    mnt  = (mnt - repmat(mean(mnt, 2), [1, nts])) ./ repmat(std(mnt, 1, 2), [1, nts]);
    sdm  = (sdm - repmat(mean(sdm, 2), [1, nts])) ./ repmat(std(sdm, 1, 2), [1, nts]);
end

CRR  = diag(mnt*sdm.')/nts;  
% the following removes any effect of polarity mismatch
mnt  = mnt .* repmat(sign(CRR), [1, nts]); CRR  = CRR.'; 
SNR  = 20 * log10(norms(sdm, 2, 2) ./ norms(mnt-sdm, 2, 2)); SNR  = SNR.';

%%
if( strcmp(dsply_opt, 'Display') )
    for i = 1 : nnd
        lgnd{i} = ['Source ' num2str(dli(i))]; 
    end
    tp = 0:1/frq:(nts-1)/frq; tp = repmat(tp.', 1, nnd);
    
    h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
    plot(tp, mnt.'), legend(lgnd), xlabel('Time (ms)'), ylabel('Moment (\muAmm)')
    
%     h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
%     for j = 1 : 3
%         subplot(3, 1, j)
%         plot(tp, mnt_cmp(:, :, j).'), legend(lgnd), xlabel('Time (ms)'), ylabel('Moment (\muAmm)'), title(['Component: ' num2str(j)])
%     end    
end

end