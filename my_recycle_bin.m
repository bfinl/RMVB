        case 'TEST' % Cam method + inflated ellipsoid
            d_disp   = 1;
            NOI_NBR  = [];
            ext_vxl_idx  = extend_ind(i_vxl, lfd_type);
            currylfd(:, ext_vxl_idx)  = mean(inverselfd(:, ext_vxl_idx, :), 3);
            
            A  = 0;  B  = 0;
            for i_nbr = 1 : length(nbr_idx)
                ext_nbr_idx  = extend_ind(nbr_idx(i_nbr), lfd_type);
                for i_mod = 1 : size(forwardlfd, 3)
                    lfd_nbr  = forwardlfd(:, ext_nbr_idx, i_mod);
                    lfd_ctr  = currylfd(:, ext_vxl_idx);
                    A  =  A + lfd_nbr.' * lfd_nbr;
                    B  =  B + lfd_nbr.' * lfd_ctr;
                end
            end
            [U, ~, V]  = svd((A^-1) * B);

            for i_nbr = 1 : length(nbr_idx)
                ext_nbr_idx  = extend_ind(nbr_idx(i_nbr), lfd_type);
                for i_mod = 1 : size(forwardlfd, 3)
                    lfd_nbr  = forwardlfd(:, ext_nbr_idx, i_mod);
                    lfd_nbr  = lfd_nbr * U * V.';
                    noi_nbr  = lfd_nbr - currylfd(:, ext_vxl_idx);
                    NOI_NBR  = [NOI_NBR noi_nbr];
                end
            end
            noisecov(:, :, i_vxl)  = find_cov(NOI_NBR, 'SC', 'do not remove mean');            
            invnoisecov  = noisecov(:, :, i_vxl) ^ -1;
            inflatefact  = max(diag(NOI_NBR.' * invnoisecov * NOI_NBR));
%             inflatefact  = chi2inv(0.975, m);
            elipsoid(:, :, i_vxl)  = noisecov(:, :, i_vxl) * inflatefact;
            [U, S, ~]    = svd(elipsoid(:, :, i_vxl));
            gammamat(:, :, i_vxl)  = U * sqrt(S) * U.';