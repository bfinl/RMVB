function [currylfd, gammamat]  = find_forward_model(LFD, TRI, type, LOC)
% this function assumes that the reference leadfield consists of the elements of the last dimension

[m, n, p]  = size(LFD);
currylfd   = mean(LFD, 3);
noisecov   = zeros(m, m , n);
gammamat   = zeros(m, m , n);
SNR   = zeros(n, 1);
NOI = [];

switch type
    case 'fxd'
        h = waitbar(0, 'Preparing to find the uncertainty...');
        for i_vxl = 1 : 1 : n
%             tri = TRI(:, find_triind(TRI, i_vxl, 1));
%             nbr = unique(tri(:));
%             nbr(nbr > n) = [];     
            
            dist = norms(LOC - repmat(LOC(:, i_vxl), [1, n]));
            nbr  = find(dist <= 10);
            nbr(nbr > n) = []; 
            
            NOI_NBR  = [];
            for i_nbr = 1 : length(nbr)            
                for i_mod = 1 : p
                    lfd_nbr  = LFD(:, nbr(i_nbr), i_mod);
                    lfd_mod  = LFD(:, i_vxl, i_mod);
                    tempppp  = pinv(lfd_nbr) * lfd_mod;
                    tempppp  = tempppp / norms(tempppp);
                    lfd_nbr  = lfd_nbr * tempppp;
                    noi_nbr  = lfd_nbr - currylfd(:, i_vxl);
                    NOI_NBR  = [NOI_NBR noi_nbr];
                end
            end
            noisecov(:, :, i_vxl)  = find_cov(NOI_NBR, 'EM', 'do not remove mean');
%             [U, S, ~]  = svd(noisecov(:, :, i_vxl) ^ -1); 
%             gammamat(:, :, i_vxl)  = (sqrt(S) * U.')^-1;
            errorstd  = sqrt(mean(diag(noisecov(:, :, i_vxl))));
            SNR(i_vxl)  = 20*log10(norms(currylfd(:, i_vxl))/sqrt(m)/errorstd);
            waitbar(i_vxl/n, h, 'Finding the uncertainty...'); 
            NOI  = [NOI NOI_NBR];
        end
        
    case 'rot'
end

% figure
% histogram(SNR,'Normalization','probability')
% title(['Ave.: ' num2str(mean(SNR)) ', STD: ' num2str(std(SNR, 1))])
% figure
% imagesc(abs(mean(noisecov, 3)))
% figure
% imagesc(abs(mean(noisecov, 3)^-1))
% figure
% for i_vxl = 1 : 10 : size(noisecov, 2)
%     imagesc(abs(noisecov(:, :, i_vxl))), colorbar
%     pause(0.030)
% end
% figure
% for i_vxl = 1 : 10 : size(noisecov, 2)
%     imagesc(abs(noisecov(:, :, i_vxl)^-1)), colorbar
%     pause(0.030)
% end

% figure 
% for i = 1 : 1
%     histogram(NOI(i, :),'Normalization','probability')
% end
% save('NOI_v4.mat', 'NOI');


IND = 1:128;
figure
hold on
for i = 1 : length(IND)

    h1 = histogram(NOI(IND(i), :),'Normalization','probability', 'FaceColor', 'r');
    Values1 = h1.Values / h1.BinWidth;
    Limits1 = linspace(h1.BinLimits(1), h1.BinLimits(2), length(Values1));
    
    plot(Limits1, Values1)
    
end


end