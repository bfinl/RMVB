function [gammamat, currylfd]  = find_uncertainty(forwardlfd, inverselfd, forwardloc, inverseloc, forwardori, inverseori, ure_opt, ext_opt, lfd_type, SRCIDXMAT)
% This function estimates the uncertainty ellipsoids for a model 
% according to section "Uncertainty Region Estimation" of the RMVM paper
% (DOI: 10.1109/TBME.2018.2859204)
% 
% Parmeters:
%     forwardlfd: a set of leadfield matrices corresponding to different uncertainty factors 
%     inverselfd: leadfield matrix for the inverse problem
%     forwardloc: 3d locations of the vertices in the forward leadfield
%     inverseloc: 3d locations of the vertices in the inverse leadfield
%     forwardori: orientation at the vertices of the forward leadfield
%     inverseori: orientation at the vertices in the inverse leadfield
%     ure_opt: uncertainty region estimation type
%     lfd_type: either 'rot' or 'fxd' for a rotational or fixed-orientation model
%     ext_opt: If "EXT", parameter deg is determined based on lfd_type (see line 28)
%     SRCIDXMAT: the indices of the points for which to estimate the uncertainty ellipsoid
% 
% Returns:
%     gammamat: the matrix to define the shape and size of the uncertainty ellipsoid (refer to the RMVB paper)
%     currylfd: the modified leadfiled to be used for solving the inverse problem
% 
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

%%
% deg: determines if the output should be a single ellipsoid for all x, y
% and z directions or there should be 3 different ellipsoids for each direction 
if(strcmp(lfd_type, 'rot')), d = 3; elseif(strcmp(lfd_type, 'fxd')), d = 1; end
if(d == 3), G = [1, 0, -1; -1, 1, 0; 0, -1, 1]; elseif(d == 1), G = 1; end
if(strcmp(ext_opt, 'EXT')), deg = d; else deg = 1; end
m  = size(inverselfd, 1); n  = size(inverselfd, 2);
currylfd  = mean(inverselfd, 3);
disp_opt  = 0;


%%
start_point1 = tic;
[~, IDX]  = find_nvoxel(forwardloc, inverseloc);
start_point2 = tic; h = waitbar(0, 'Preparing to find the uncertainty...');
parfor i_vxl = 1 : 1 : n/d
% for i_i_vxl = 1 : length(SRCIDXMAT)
%     i_vxl  = SRCIDXMAT(i_i_vxl);
    
    %=====================================================================%%=====================================================================%
    nbr_idx  = find(IDX == i_vxl);
    if(isempty(nbr_idx)) % this can happen, since there could be repititious locations
        [~, i_vxl_peer]  = find_nvoxel(inverseloc(:, i_vxl), inverseloc);
        nbr_idx  = find(IDX == i_vxl_peer); 
        if(isempty(nbr_idx))
            loc  = find_nvoxel(inverseloc(:, i_vxl), forwardloc);
            r    = norm(loc-inverseloc(:, i_vxl)) + 1;
            dist = norms(forwardloc - repmat(inverseloc(:, i_vxl), [1, size(forwardloc, 2)]));
            nbr_idx  = find(dist <= r);
        end
    end
    
    %=====================================================================%%=====================================================================%
    if(strcmp(lfd_type, 'fxd'))
        tmpori   = find_angle(inverseori(:, i_vxl), forwardori(:, nbr_idx), '[0, 180]');
        anglestd(i_vxl, 1)  = std(tmpori, 1);
        nbr_idx  = nbr_idx(tmpori < 180);
    end
            
    %=====================================================================%%=====================================================================%
    lfd_ctr  = 0; NOI_NBR  = []; LFD_NBR  = [];
    ext_vxl_idx  = extend_ind(i_vxl, lfd_type);
    for i_nbr = 1 : length(nbr_idx)
        ext_nbr_idx  = extend_ind(nbr_idx(i_nbr), lfd_type);
        for i_mod = 1 : size(forwardlfd, 3)
            lfd_nbr  = forwardlfd(:, ext_nbr_idx, i_mod);
            switch ure_opt % uncertainty region estimation type
                case 'EMP'
                    inflation_opt = 'incallneighbours';
                    lfd_ctr  = currylfd(:, ext_vxl_idx);
                
                case 'EMP-v2'
                    inflation_opt = 'incallneighbours';
                    lfd_ctr = lfd_ctr + lfd_nbr/length(nbr_idx)/size(forwardlfd, 3);

                case 'CAM'
                    inflation_opt = 'incallneighbours';
                    lfd_ctr  = currylfd(:, ext_vxl_idx);
                    lfd_mod  = inverselfd(:, ext_vxl_idx, i_mod);
                    F  = pinv(lfd_nbr) * lfd_mod; F  = F / norm(F, 'fro');
                    lfd_nbr  = lfd_nbr * F;

                case 'AMIR'
                    inflation_opt = 'confidencemargin';
                    lfd_ctr  = currylfd(:, ext_vxl_idx);
                    lfd_mod  = inverselfd(:, ext_vxl_idx, i_mod);
                    [U, ~, V]  = svd(pinv(lfd_nbr) * lfd_mod); F = U * V.';
                    lfd_nbr  = lfd_nbr * F;

                case 'AMIR-v2'
                    inflation_opt = 'confidencemargin';
                    lfd_ctr  = currylfd(:, ext_vxl_idx);
                    [Ud, ~, ~] = svds(lfd_ctr * G, d); 
                    [U, ~, V]  = svd(pinv(Ud.'*lfd_nbr) * (Ud.'*lfd_ctr)); F = U * V.';
                    lfd_nbr  = lfd_nbr * F;
                    
                case 'TEST'
                    inflation_opt = 'incallneighbours';
                    lfd_ctr = lfd_ctr + lfd_nbr/length(nbr_idx)/size(forwardlfd, 3);

            end 
            LFD_NBR  = [LFD_NBR lfd_nbr];
        end
    end
    [A, B, alpha]  = estimate_ellipsoid_v2(lfd_ctr, LFD_NBR, deg, d, inflation_opt, disp_opt);
    elipsoid_cell{1, i_vxl}  = A;    gammamat_cell{1, i_vxl}  = B;   Alpha{1, i_vxl}  = alpha;
    
    %=====================================================================%%=====================================================================%
%     percentage  = i_vxl*d/n;
%     waitbar(percentage, h, ['Finding the uncertainty, ' datestr((1-percentage)/percentage*toc(start_point2)/60/60/24, 'HH:MM:SS') ' left...']); 
        
end
toc(start_point1)


%%
Alpha     = cat(1, Alpha{:});
elipsoid  = cat(3, elipsoid_cell{:});
gammamat  = cat(3, gammamat_cell{:});
[elipsoid, gammamat]  = rescale_ellipsoid(elipsoid, gammamat, deg, d, Alpha, 'mean', currylfd, anglestd);


%%
% figure, imagesc(abs(mean(elipsoid, 3))), colorbar, title(ure_opt)
% find_spatial_orientation(abs(mean(elipsoid, 3))) 

            
end


function [elipsoid, gammamat, alpha]  = estimate_ellipsoid_v2(lfd_ctr, LFD_NBR, deg, d, inflation_option, disp_option)
%=========================================================================%%=========================================================================%
alpha     = zeros(1, deg);
elipsoid  = zeros(size(LFD_NBR, 1), size(LFD_NBR, 1), deg);
gammamat  = zeros(size(LFD_NBR, 1), size(LFD_NBR, 1), deg);
%=========================================================================%%=========================================================================%
NOI_NBR  = LFD_NBR - repmat(lfd_ctr, [1, size(LFD_NBR, 2)/size(lfd_ctr, 2)]);
for i_cmp = 1 : deg
    NOI_NBR_CMP  = NOI_NBR(:, i_cmp:deg:end);
    noisecov  = find_cov(NOI_NBR_CMP, 'SC', 'do not remove mean');
    switch inflation_option
        case 'covariancematrix'
            inflatefact  = 1;
            INDINF  = [];
            
        case 'incallneighbours' % includes all neighbour points
            [inflatefact, indinf]  = max(diag(NOI_NBR_CMP.' * (noisecov ^ -1) * NOI_NBR_CMP)); 
            INDINF(i_cmp) = 3*(indinf-1) + i_cmp;
            
        case 'confidencemargin' % based on confidence margin
            inflatefact  = chi2inv(0.975, size(NOI_NBR, 1));
            INDINF  = [];
            
        otherwise
            display('Wrong inflation option!!!')
    end
    elipsoid(:, :, i_cmp)  = noisecov * inflatefact;
    [U, S, ~]              = svd(elipsoid(:, :, i_cmp));  
    gammamat(:, :, i_cmp)  = U * (S.^0.5) * U.';       
end
elipsoid  = repmat(elipsoid, [1, 1, d/size(elipsoid, 3)]);
gammamat  = repmat(gammamat, [1, 1, d/size(gammamat, 3)]);
%=========================================================================%%=========================================================================%
for i_cmp = 1 : d
    alpha(i_cmp)  = ...
        1 ./ sqrt(lfd_ctr(:, i_cmp).' * (elipsoid(:, :, i_cmp)^-1) * lfd_ctr(:, i_cmp));
end
%=========================================================================%%=========================================================================%
if(disp_option)
    my_minivisualization(lfd_ctr, gammamat, d, 'Display', LFD_NBR, INDINF);
end
end


function [elipsoid, gammamat] = rescale_ellipsoid(elipsoid, gammamat, deg, d, Alpha, rescale_opt, currylfd, anglestd)
%=========================================================================%
rAlpha = Alpha(Alpha < 1);
switch rescale_opt
    case 'mean'
        scale_thresh  = mean(rAlpha);
    case 'median'
        scale_thresh  = median(rAlpha);
    case '90-percentile'
        scale_thresh  = prctile(rAlpha, 90);
end
%=========================================================================%
Alpha_new  = Alpha;
for i = 1 : d/deg : length(Alpha)
    Alpha_new(i:i+d/deg-1)  = max(Alpha(i:i+d/deg-1));
end
%=========================================================================%
parfor i = 1 : size(gammamat, 3)
    if(Alpha_new(i) > scale_thresh)
        scale_factor      = scale_thresh / Alpha_new(i);
        gammamat(:, :, i) = gammamat(:, :, i) * scale_factor;
        elipsoid(:, :, i) = elipsoid(:, :, i) * (scale_factor^2);
    end
end
fprintf('Scaling Threshold is %f\n', scale_thresh);
%=========================================================================%
h = figure;  h = ApplyProperties(h, 'Customized-v1');
plot(Alpha, anglestd, 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.25, 'MarkerFaceColor', 'g'), 
xlabel('\beta'), ylabel('Angle STDev (Deg)')
% figure, scatter(Alpha, norms(currylfd)), xlabel('\alpha'), ylabel('Leadfield Norm')
myper  = [10:10:100 25, 75];

h = figure;  h = ApplyProperties(h, 'Customized-v1');
plot(anglestd, Alpha, 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.25, 'MarkerFaceColor', 'g'), 
xlabel('Angle STDev (Deg)'), ylabel('\beta')
% figure, scatter(Alpha, norms(currylfd)), xlabel('\alpha'), ylabel('Leadfield Norm')
myper  = [10:10:100 25, 75];

h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
histogram(Alpha, 100, 'EdgeAlpha', 0, 'FaceColor', [20, 20, 20]/256); lgnd{1} = 'Hist'; hold on, xlabel('\beta'), ylabel('#')
line([mean(rAlpha), mean(rAlpha)], [0, 1200], 'LineWidth', 2, 'Color', [1 0.2 0]); lgnd{2} = ['Mean: ' num2str(mean(rAlpha))]; 
% line([median(rAlpha), median(rAlpha)], [0, 1200], 'LineWidth', 2, 'Color', [0.25 0.5 0]); lgnd{3} = ['Median: ' num2str(median(rAlpha))]; 
legend(lgnd)

h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
histogram(Alpha, 100, 'Normalization', 'pdf'); lgnd{1} = 'Hist'; hold on
stem(mean(rAlpha), 1.457, ':diamond', 'LineWidth', 2); lgnd{2} = ['mean: ' num2str(mean(rAlpha))]; 
for i = 1 : length(myper)
    p = prctile(rAlpha, myper(i)); stem(p, 1, 'LineWidth', 2);
    lgnd{i+2} = ['percentile %' num2str(myper(i)) ': ' num2str(p)]; 
end
legend(lgnd)

h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
h = histogram(Alpha, 100, 'Normalization', 'cdf'); lgnd{1} = 'Hist'; hold on
stem(mean(rAlpha), 1, ':diamond', 'LineWidth', 2); lgnd{2} = ['mean: ' num2str(mean(rAlpha))]; 
for i = 1 : length(myper)
    p = prctile(rAlpha, myper(i)); stem(p, 1, 'LineWidth', 2);
    lgnd{i+2} = ['percentile %' num2str(myper(i)) ': ' num2str(p)]; 
end
legend(lgnd)


end


% function [elipsoid, gammamat]  = estimate_ellipsoid_v1(lfd_ctr, LFD_NBR, deg, inflation_option, disp_option, varargin)
% alpha     = 0.01;
% elipsoid  = zeros(size(LFD_NBR, 1), size(LFD_NBR, 1), deg);
% gammamat  = zeros(size(LFD_NBR, 1), size(LFD_NBR, 1), deg);
% [Unorm, ~, ~]  = svds(lfd_ctr, size(lfd_ctr, 2)); Uproj  = null(Unorm.'); 
% NOI_NBR  = LFD_NBR - repmat(lfd_ctr, [1, size(LFD_NBR, 2)/size(lfd_ctr, 2)]);
% for i_cmp = 1 : deg
%     NOI_NBR_CMP  = NOI_NBR(:, i_cmp:deg:end);
%     noisecov  = find_cov(NOI_NBR_CMP, 'SC', 'do not remove mean');
%     switch inflation_option
%         case 'covariancematrix'
%             inflatefact  = 1;
%             INDINF  = [];
%             
%         case 'incallneighbours' % includes all neighbour points
%             [inflatefact, indinf]  = max(diag(NOI_NBR_CMP.' * (noisecov ^ -1) * NOI_NBR_CMP)); 
%             INDINF(i_cmp) = 3*(indinf-1) + i_cmp;
%             
%         case 'confidencemargin' % based on confidence margin
%             inflatefact  = chi2inv(0.975, size(NOI_NBR, 1));
%             INDINF  = [];
%             
%         otherwise
%             display('Wrong inflation option!!!')
%     end
%     elipsoid(:, :, i_cmp)  = noisecov * inflatefact;
%     [U, S, ~]              = svd(Uproj * (Uproj.') * elipsoid(:, :, i_cmp) * Uproj * (Uproj.'));  
%     S                      = (1-alpha) * S + alpha * mean(diag(S)) * eye(size(S));
%     gammamat(:, :, i_cmp)  = U * (S.^0.5) * U.';   
% end
% if(disp_option)
%     d  = varargin{1};
%     my_minivisualization(lfd_ctr, repmat(gammamat, [1, 1, d/size(gammamat, 3)]), d, 'Display', LFD_NBR, INDINF);
% end
% end






