% ratiolabel = {'0.0400', '0.0500', '0.0667'};
% tic
% ind  = cell(2, length(ratiolabel));
% for i = 2 : length(ratiolabel)-1
%     load(['Resource\lfd_BEM-' ratiolabel{i} '_CTX-1.0mm-fxd_EEG-Biosemi-128.mat'])
%     lfd0 = currylfd; loc0 = curryloc; 
%     load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128.mat')
%     lfd1 = currylfd; loc1 = curryloc; 
% 
%     [~, ind{1, i}] = find_nvoxel(loc0, loc1);
% %     [~, ind{2, i}] = find_nvoxel(lfd0, lfd1);
% end
% % save('matlab_estimate_snr.mat')
% toc

%%
ratiolabel = {'0.0400', '0.0500', '0.0667'};
load('matlab_estimate_snr.mat')
for IND = 1 : 1
    NOI = [];
    for i = 1 : length(ratiolabel)
        load(['Resource\lfd_BEM-' ratiolabel{i} '_CTX-1.0mm-fxd_EEG-Biosemi-128.mat'])
        lfd0 = currylfd; loc0 = curryloc; 
        load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128.mat')
        lfd1 = currylfd; loc1 = curryloc; 

        lfd1 = lfd1(:, ind{IND, i}); 
        loc1 = loc1(:, ind{IND, i});

        noi = lfd1-lfd0;
        dis = norms(loc0-loc1);
        snr = 20*log10(norms(lfd0) ./ norms(noi));
        figure, histogram(snr,'Normalization','probability')
        title(['Ratio: ' ratiolabel{i} ', Ave.: ' num2str(mean(snr)) ', STD: ' num2str(std(snr, 1))])
        figure, histogram(dis,'Normalization','probability'), xlabel('Distance (mm)')
        figure, imagesc(abs(noi * noi.' / size(noi, 2)))
        figure, hold on
        NOI = [NOI, noi];
    end
    for i = 1 : 1
            histogram(NOI(i, :),'Normalization','probability', 'FaceColor', 'r')
    end
end
save('NOI_REF.mat', 'NOI');






%%
% ratiolabel = {'0.0400', '0.0500', '0.0667'};
% tic
% ind  = cell(2, length(ratiolabel));
% for i = 2 : length(ratiolabel)-1
%     load(['Resource\lfd_BEM-' ratiolabel{i} '_CTX-3.0mm-fxd_EEG-Biosemi-128.mat'])
%     lfd0 = currylfd; loc0 = curryloc; ori0 = curryori;
% %     load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128.mat')
%     load(['Resource\lfd_BEM-' ratiolabel{i} '_CTX-5.0mm-fxd_EEG-Biosemi-128.mat'])
%     lfd1 = currylfd; loc1 = curryloc; ori1 = curryori; 
% 
%     [~, ind] = find_nvoxel(loc1, loc0);
%     figure
%     x = ori1; y = ori0(:, ind);
%     h = histogram(2 * atan(norms(x - y) ./ norms(x + y)) * 180 / pi, 'Normalization','probability')
%     Values = h.Values / h.BinWidth;
%     Limits = linspace(h.BinLimits(1), h.BinLimits(2), length(Values));
%     
%     plot(Limits, Values)
%     
% end
% toc

load('Resource\lfd_BEM-[0.0400-0.0667-11]_CTX-5.0mm-fxd_EEG-Biosemi-128.mat')
loc1 = curryloc; 

r = [];
for i = 1 : size(loc1, 2)
    loc = draft(loc1(:, i), loc1);
    r = [r  norm(loc-loc1(:, i))];
end
figure
histogram(r, 800)


