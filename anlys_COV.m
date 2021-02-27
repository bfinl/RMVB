function [rel_err]  = anlys_COV(tru_epl, est_epl, tru_cov, currycdr, curryloc)

rel_err  = 0;
tru_epi  = ppatches(find_nvoxel_interface(tru_epl, curryloc), 'row');
est_epi  = ppatches(find_nvoxel_interface(est_epl, curryloc), 'row');
cur_epi  = [tru_epi setdiff(est_epi, tru_epi)];

ref_cov  = zeros(length(cur_epi));
ref_cov(1:length(tru_epi), 1:length(tru_epi)) = tru_cov;
est_cdr  = revise_cdr_ori(currycdr(:, cur_epi, :), inf);
est_cdr  = squeeze(est_cdr(4, :, :));
est_cov  = find_cov(est_cdr, 'EM');


figure, colormap hot
ref_cov  = ref_cov/norm(ref_cov, 'fro');
est_cov  = est_cov/norm(est_cov, 'fro');
fprintf('Covariance Recovery SNR : %f\n', 20*log10(norm(est_cov-ref_cov, 'fro')/norm(ref_cov, 'fro')))
subplot(1, 2, 2), imagesc(abs(est_cov)/max(max(abs(est_cov)))), colorbar, title('Estimated Covariance')
subplot(1, 2, 1), imagesc(abs(ref_cov)/max(max(abs(ref_cov)))), colorbar, title('Ground Truth Covariance')




end


function epi  = find_nvoxel_interface(epl, curryloc)

epi = cell(1, length(epl));
for i = 1 : length(epl)
    [~, epi{i}]  = find_nvoxel(epl{i}, curryloc);
end

end