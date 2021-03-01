% This program is a free software for academic research: you can redistribute 
% it and/or modify it for non-commercial uses, under the license terms provided 
% with the package at the GitHub page where this package is downloaded. 
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the License for more details. This program is 
% for research purposes only. This program CAN NOT be used for commercial purposes. 
% This program SHOULD NOT be used for medical purposes. The authors WILL NOT 
% be responsible for using the program in medical conditions.
% ==========================================

function [CDR, Gamma, C_bl]  = solve_inverse_problem(Phi_bl, Phi_az, LFD, ULFD, ORI, TRI, ANI, lfd_type, cdr_type, nnd, varargin)
% This function solves the inverse problem to generate the underlying source activity
% corresponding to a given set of electrode recordings using different reconstruction techniques. 
% 
% Parmeters:
%     Phi_bl: the baseline electrode recordings
%     Phi_az: the electrode recordings when the sources are active (active zone)
%     LFD: the inverse leadfield matrix
%     ORI: the orientation of source dipoles at each 3D location of the brain. 
%     TRI: the traingular mesh of the brain.
%     lfd_type: "fxd" or "rot" for a fixed or rotational orientation modeling of the sources
%     nnd: number of active nodes in the source
% 
% Returns:
%     CDR: The reconstructed source activity over 3D locations and time
%         
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017

Gamma = [];
C_bl  = [];
ROIA  = 0;

switch lfd_type 

%=========================================================================%%=========================================================================%
    case 'rot'
        
        switch cdr_type

            %=============================================================%%=============================================================%
            case 'LS'
 
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, S, V]  = svd(A, 'econ');
                J  = V * (S^-1) * U.' * b;
                CDR  = currize(J, [], 'rot');


            %=============================================================%%=============================================================%
            case 'LS-ROI'
 
                C_bl  = find_cov(Phi_bl, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', zeros(1, nnd));
                ExtIDX   = extend_ind(IDX, 'rot');
                [A, b]  = mywhitening(LFD(:, ExtIDX), Phi_az, C_bl);
                [U, S, V]  = svd(A, 'econ');
                J  = zeros(size(LFD, 2), size(b, 2));
                J(ExtIDX, :)  = V * (S^-1) * U.' * b;
                CDR  = currize(J, [], 'rot');
                
                
            %=============================================================%%=============================================================%
            case 'MN-DCP'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta);  
                lambda  = lambda ^ 2;
                T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
                J  = T * b;
                CDR  = currize(J, [], 'rot');            

                
            %=============================================================%%=============================================================%
            case 'WMN-DCP'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                w  = myblockaverage((norms(A).^2).', 3).^-1; 
                A  = A * mysparseblockdiag(w.^0.5);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta);  
                lambda  = lambda ^ 2;
                T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
                J  = mysparseblockdiag(w.^0.5) * (T * b); 
                CDR  = currize(J, [], 'rot'); 

                
            %=============================================================%%=============================================================%
            case 'sLORETA-DCP'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta);  
                lambda  = lambda ^ 2;
                T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
                C  = zeros(3, 3, size(T, 1)/3);
                for i = 1 : 3 : size(T, 1)
                    C(:, :, floor(i/3)+1)  = T(i:i+2, :) * A(:, i:i+2);
                end
                iC  = my3by3inv(C);
                J  = T * b;
                J  = reshape(J, [1, 3, size(J, 1)/3, size(J, 2)]);
                CDR(4, :, :)  = squeeze(mmat(mmat(J, iC), permute(J, [2, 1, 3, 4])));
                CDR(1:3, :, :)  = squeeze(J) ./ repmat(norms(squeeze(J), 2, 1), [3, 1, 1]);
                
                
            %=============================================================%%=============================================================%
            case 'sLORETA-DCP-TC'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                iC  = zeros(3, 3, size(A, 2)/3, size(b, 2));
                C   = zeros(3, 3, size(A, 2)/3);
                J   = zeros(size(A, 2), size(b, 2)); 
                for t = 1 : size(b, 2)
                    [~, lambda] = discrep(U, s, V, b(:, t), delta);  
                    lambda  = lambda ^ 2;
                    T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
                    for i = 1 : 3 : size(T, 1)
                        C(:, :, floor(i/3)+1)  = T(i:i+2, :) * A(:, i:i+2);
                    end
                    iC(:, :, :, t)  = my3by3inv(C);
                    J(:, t)  = T * b(:, t);
                end 
                J  = reshape(J, [1, 3, size(J, 1)/3, size(J, 2)]);
                CDR(4, :, :)  = squeeze(mmat(mmat(J, iC), permute(J, [2, 1, 3, 4])));
                CDR(1:3, :, :)  = squeeze(J) ./ repmat(norms(squeeze(J), 2, 1), [3, 1, 1]);
                
                
            %=============================================================%%=============================================================%
            case 'MUSIC'
                
                C_az  = find_cov(Phi_az, 'SC');
                [U_s, ~, ~]  = svds(C_az, nnd);                
                CDR  = zeros(4, size(LFD, 2)/3, size(Phi_az, 2));
                count  = 0;
                for i = 1 : 3 : size(LFD, 2)
                    count  = count + 1;
                    CDR(4, count, :)  = max(eig((ULFD(:, i:i+2).') * U_s * (U_s.') * ULFD(:, i:i+2)));
                end
                
                
            %=============================================================%%=============================================================%
            case 'LCMV'

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : 3 : size(LFD, 2)
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                end

                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND'

                iC_bl  = find_cov(Phi_bl, 'SC')^-1;
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : 3 : size(LFD, 2)
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                    T(i:i+2, :) = T(i:i+2, :) / sqrt(trace((LFD(:, i:i+2).' * iC_bl * LFD(:, i:i+2)) ^ -1));
                end

                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-iMAP'

                iC_bl  = find_cov(Phi_bl, 'SC')^-1;
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                varN   = zeros(size(LFD, 2)/3, 1);
                count  = 0;
                for i = 1 : 3 : size(LFD, 2)
                    count  = count + 1;
                    varN(count) = trace((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1); 
                    varN(count) = varN(count) / trace((LFD(:, i:i+2).' * iC_bl * LFD(:, i:i+2)) ^ -1);
                end
                
                CDR  = zeros(4, size(LFD, 2)/3, size(Phi_az, 2));
                CDR(4, :, : )  = repmat(varN, [1, 1, size(Phi_az, 2)]);

                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-DN'
 
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                iC_bl  = C_bl^-1;
                iC_az  = C_az^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : 3 : size(LFD, 2)
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                    T(i:i+2, :) = T(i:i+2, :) / sqrt(trace((LFD(:, i:i+2).' * iC_bl * LFD(:, i:i+2)) ^ -1));
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                                
                
            %=============================================================%%=============================================================%
            case 'Champagne-D'

                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = ones(size(LFD, 2), 1);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                epsil = 1e-5; iter_max = 100; ofval = 0;
                for iter = 1 : iter_max                     
                    P  = LFD.' * iS_az;
                    X  = P * C_sr;
                    Z  = dot(P, LFD.', 2);
                    gamma  = gamma .* norms(X, 2, 2) .* (Z.^-0.5); 
                    Gamma  = mysparseblockdiag(gamma);
                    S_az   = LFD * Gamma * LFD.' + C_bl;
                    iS_az  = S_az ^ -1;
                    fval  = trace(C_az * iS_az) + sum(log(eig(S_az)));
                    if(abs((fval - ofval)/fval) <= epsil), break, end
                    ofval  = fval;
                end   
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'Champagne-S' % the same as SBL

                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = ones(size(LFD, 2), 1);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                epsil = 1e-5; iter_max = 100; ofval = 0;
                for iter = 1 : iter_max     
                    iter
                    P  = LFD.' * iS_az;
                    X  = P * C_sr;
                    Z  = dot(P, LFD.', 2);
                    gamma  = gamma .* (myblockaverage(norms(X, 2, 2).^2, 3).^0.5) .* (myblockaverage(Z, 3).^-0.5);
                    Gamma  = mysparseblockdiag(gamma);
                    S_az   = LFD * Gamma * LFD.' + C_bl;
                    iS_az  = S_az ^ -1;
                    fval  = trace(C_az * iS_az) + sum(log(eig(S_az)));
                    if(abs((fval - ofval)/fval) <= epsil), break, end
                    ofval  = fval;
                end   
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 

                
            %=============================================================%%=============================================================%
            case 'Champagne-M'

                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = repmat(eye(3, 3), [size(LFD, 2)/3, 1]) + 1e-5 * randn(size(LFD, 2), 3);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                epsil = 1e-5; iter_max = 100; ofval = 0;
                for iter = 1 : iter_max 
                    P  = LFD.' * iS_az;
                    X  = Gamma * (P * C_sr);
                    for i = 1 : 3 : size(LFD, 2)
                        Z  = P(i:i+2, :) * LFD(:, i:i+2); [hZ, ihZ]  = myx2factrzn(Z);
                        gamma(i:i+2, :)  = ihZ * myx2factrzn(hZ * X(i:i+2, :) * X(i:i+2, :).' * hZ) * ihZ;
                    end
                    Gamma  = mysparseblockdiag(gamma);
                    S_az   = LFD * Gamma * LFD.' + C_bl;
                    iS_az  = S_az ^ -1;
                    fval   = trace(C_az * iS_az) + sum(log(eig(S_az)));
                    if(abs((fval - ofval)/fval) <= epsil), break, end
%                     abs((fval - ofval)/fval)
%                     fval-ofval
                    ofval  = fval;
                    costt(iter) = fval;
                end   
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');  
                
                figure, plot(costt);


            %=============================================================%%=============================================================%
            case 'SBL-IR-L2'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                w  = ones(size(LFD, 2), 1);
                
                epsil = 1e-6; iter_max = 100; alpha = 1e-3;
                for iter = 1 : iter_max
                    iter 
                    [J, lambda] = mywmn(A, b, w.^-1);
                    alpha = alpha * 0.9;
                    W   = mysparseblockdiag(w.^-1);
                    Z   = dot(A.'*(A*W*A.' + alpha*eye(size(A, 1)))^-1, A.', 2);
                    nw  = (norms(J, 2, 2).^2 + w.^-1 - (w.^-2).*Z).^-1; 
                    if(max(abs(nw-w)) < epsil), break, end
                    w   = nw;
                end                                   
                
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'GLAII'
                
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [~, azp_ind]  = max(var(b));
                [A, L]  = revise_grp(A, ANI, 'rot', 'Forward');

%                 w  = myblockaverage((norms(A).^2).', 3).^-1; 
%                 A  = A * mysparseblockdiag(w.^0.5);
                [z, ~]  = group_lasso(A, b(:, azp_ind), 1, L, 10, 1.4);
%                 z  = mysparseblockdiag(w.^0.5) * z; 

                [J, ~]  = revise_grp(z, ANI, 'rot', 'Inverse');
                CDR  = currize(J, [], 'rot', size(Phi_az, 2)); 
               
              
            %=============================================================%%=============================================================%
            case 'GLAII-CVX'
                
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                delta   = sqrt(chi2inv(0.975, size(b, 1)));
                [~, azp_ind]  = max(var(b));
                J  = glasso_cvx(A, b(:, azp_ind), delta, ANI, 'regular', [], 'rot');
%                 J  = glasso_cvx(A, b(:, azp_ind), delta, ANI, 'sparse GL', 0, 'rot');
                CDR  = currize(J, [], 'rot', size(Phi_az, 2));
                
              
            %=============================================================%%=============================================================%
            case 'ENet'
                
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [~, azp_ind]  = max(var(b));
%                 [J, ~] = elastic_net_lrn(A, b(:, azp_ind), 1, 1);
%                 J  = stepwisefit(A, b(:, azp_ind));
                [J inf]  = lassoglm(A, b(:, azp_ind), 'normal', 'Alpha', 0.9);
                
%                 lassoPlot(J,FitInfo,'plottype','CV');
                
%                 J  = repmat(J, [1, size(Phi_az, 2)]);
                CDR  = currize(J, [], 'rot');
                
                
            %=============================================================%%=============================================================%
            case 'BAYES'
                
                source  = varargin{1};
                LOC  = varargin{2};
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                Gamma  = sparse(size(LFD, 2), size(LFD, 2));
                tru_epi  = ppatches(find_nvoxel_interface(source.epl, LOC), 'row');
                cov = zeros(3*length(tru_epi));
                for i = 1 : length(tru_epi)
                    for j = 1 : length(tru_epi)
                        cov(3*i-2:3*i, 3*j-2:3*j)  = ORI(:, tru_epi(i))*ORI(:, tru_epi(j)).'*source.epc(i, j);
                    end
                end
                tru_epi  = extend_ind(tru_epi, 'rot');
                Gamma(tru_epi, tru_epi)  = cov;
                T  = Gamma * A.' * (A * Gamma * A.' + eye(size(b, 1))) ^ -1;
                J  = T * b; 
                CDR  = currize(J, [], 'rot');

                
            %=============================================================%%=============================================================%
            case 'ChM-ANI' % extent sbl

                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = ones(size(LFD, 2), 1);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                epsil = 1e-5; iter_max = 100; ofval = 0;
                for iter = 1 : iter_max 
                    P  = LFD.' * iS_az;
                    X  = Gamma * (P * C_sr);
                    for i  = 1 : length(ANI)
                        I  = extend_ind(double(ANI{i}), 'rot');
                        Z  = P(I, :) * LFD(:, I); [hZ, ihZ]  = myx2factrzn(Z);
                        Gamma(I, I)  = ihZ * myx2factrzn(hZ * (X(I, :)*X(I, :).') * hZ) * ihZ;
                    end
                    S_az  = LFD * Gamma * LFD.' + C_bl;
                    iS_az  = S_az ^ -1;
                    fval  = trace(C_az * iS_az) + sum(log(eig(S_az)));
                    if(abs((fval - ofval)/fval) <= epsil), break, end
                    ofval  = fval;
                end   
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');
                
                
            %=============================================================%%=============================================================%
            case 'SBL-DDP' % SBL with data driven parcellation
                
                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = ones(size(LFD, 2), 1);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                [SDI, TAG]  = generate_seed_points(Phi_az, LFD, [], TRI, 4, 512, 'GROVA', 'rot');
                [DDP, ~]  = generate_parcellation(SDI, TAG, TRI, size(LFD, 2)/3);

                for itre_ddp  = 1 : 1             

                        epsil = 1e-5; iter_max = 20; ofval = 0;
                        for iter = 1 : iter_max 
                            P  = LFD.' * iS_az;
                            X  = P * C_sr;
                            Z  = dot(P, LFD.', 2);
                            gamma  = gamma .* (mygblockaverage(norms(X, 2, 2).^2, DDP).^0.5) .* (mygblockaverage(Z, DDP).^-0.5);
                            Gamma  = mysparseblockdiag(gamma);
                            S_az   = LFD * Gamma * LFD.' + C_bl;
                            iS_az  = S_az ^ -1;
                            fval  = trace(C_az * iS_az) + sum(log(eig(S_az)));
                            if(abs((fval - ofval)/fval) <= epsil), break, end
                            ofval  = fval;
                            costt(iter) = fval;
                        end   

                        T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                        J  = T * Phi_az;
                        CDR  = currize(J, [], 'rot'); 
                        
                        figure, plot(costt);
                        
                        MSP  = myblockaverage(gamma, 3); MSP  = MSP(1:3:end);
%                         MSP  = norms(squeeze(CDR(4, :, :)), 2, 2).^2;
%                         [MSP, ~, ~]  = svds(squeeze(CDR(4, :, :)), nnd);
                        [SDI, TAG]  = generate_seed_points(Phi_az, LFD, [], TRI, 4, 512, 'EXTNL', 'rot', MSP);
                        [DDP, ~]  = generate_parcellation(SDI, TAG, TRI, size(LFD, 2)/3);
                        
                end
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
%                 J  = repmat(gamma, [1, size(Phi_az, 2)]);
                CDR  = currize(J, [], 'rot');
                
                
            %=============================================================%%=============================================================%
            case 'SBL-ANI' % SBL with anatomical information incorporated

                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                gamma  = ones(size(LFD, 2), 1);
                Gamma  = mysparseblockdiag(gamma);
                S_az   = LFD * Gamma * LFD.' + C_bl;                
                iS_az  = S_az ^ -1;
                
                epsil = 1e-5; iter_max = 100; ofval = 0;
                for iter = 1 : iter_max 
                    iter
                    P  = LFD.' * iS_az;
                    X  = Gamma * P * C_sr;
                    Z  = dot(P, LFD.', 2);
                    gamma  = (mygblockaverage(norms(X, 2, 2).^2, ANI).^0.5) .* (mygblockaverage(Z, ANI).^-0.5);
                    Gamma  = mysparseblockdiag(gamma);
                    S_az   = LFD * Gamma * LFD.' + C_bl;
                    iS_az  = S_az ^ -1;
                    fval  = trace(C_az * iS_az) + sum(log(eig(S_az)));
                    if(abs((fval - ofval)/fval) <= epsil), break, end
                    ofval  = fval;
                end   
                
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
%                 J  = repmat(gamma, [1, size(Phi_az, 2)]);
                CDR  = currize(J, [], 'rot');

                
            %=============================================================%%=============================================================%
            case 'SBL-SVD' % SBL in source and variation domain
                
                NHD   = cell(size(LFD, 2)/3, 1);
                for i = 1 : size(LFD, 2)/3
                    ind  = TRI(:, find_triind(TRI, i, 1)); ind  = unique(ind(:)); 
                    ind(ind == i | (ind > size(LFD, 2)/3))  = [];  NHD{i}  = ind;
                end
                
%                 Alpha  = linspace(0.05, 0.015, 11);
                Alpha  = 0.001;
%                 Lambda = logspace(log10(1e-5), log10(1e5), 200);
                Lambda = 1e7;
                
                C_bl   = find_cov(Phi_bl, 'SC');
                C_az   = find_cov(Phi_az, 'SC');
                C_sr   = myxxtfactrzn(C_az);
                swei   = ones(size(LFD, 2), 1); % source weight
                vwei   = ones(size(LFD, 2), 1); % variation weight
                
                n_w  = 1; n_a  = length(Alpha); n_l  = length(Lambda);
                h = waitbar(0, 'Tuning regularization parameters...');
                for i_w = 1 : n_w;
                    for i_a = 1 : length(Alpha)
                        for i_l = 1 : length(Lambda)
                            
                            gamma  = ones(size(LFD, 2), 1);
                            Gamma  = mysparseblockdiag(gamma);
                            S_az   = LFD * Gamma * LFD.' + C_bl; 
                            iS_az  = S_az ^ -1;

                            alpha = Alpha(i_a), lambda = Lambda(i_l),
                            epsil = 1e-5; iter_max = 100; ofval = 0;
                            for iter = 1 : iter_max
                                iter
                                P  = LFD.' * iS_az;
                                X  = Gamma * (P*C_sr);
                                Z  = dot(P, LFD.', 2);
                                b  = (lambda * (1-alpha) * sqrt(3) / 2) * vwei;
                                x  = (3*myblockaverage(norms(X, 2, 2).^2, 3)).^0.5;
                                a  = 3*myblockaverage(Z, 3) + lambda*alpha*sqrt(3)*swei;
                                gamma  = solve_sparse_problem(gamma(1:3:end), x(1:3:end), a(1:3:end), b(1:3:end), NHD, 3);
                                Gamma  = mysparseblockdiag(gamma);
                                S_az   = LFD * Gamma * LFD.' + C_bl;
                                iS_az  = S_az ^ -1;

                                [fval, val1, val2]  = find_fval(NHD, alpha, lambda, 'REG', gamma(1:3:end), C_az, S_az, iS_az);
%                                 [fval, val1, val2]  = find_fval(NHD, alpha, lambda, 'CVX', gamma, X, Z, 3);
                                if(abs((fval - ofval)/fval) <= epsil), break, end
                                ofval  = fval;
                                costt(i_a, i_l) = fval; costx(i_a, i_l) = val1; costy(i_a, i_l) = val2;
                            end  
                            waitbar(((i_w-1)*n_a*n_l + (i_a-1)*n_l + i_l) / (n_w*n_a*n_l), h, 'Tuning regularization parameters...'); 
                            
                        end
                    end
                    [swei, vwei]  = update_weights(gamma(1:3:end), NHD, 3);
                end
                close(h)

                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');  
                
                myregplot(costx, costy, Alpha, Lambda);
%                 save('regparam.mat', 'costt', 'costx', 'costy', 'Alpha', 'Lambda')


            %=============================================================%%=============================================================%
            case 'MyCOV' % extent sbl     
                TH  = 0.01;
                [~, Gamma, C_bl]  = solve_inverse_problem(Phi_bl, Phi_az, LFD, ULFD, ORI, TRI, ANI, 'rot', 'SBL-ANI', nnd, varargin);
                gamma  = diag(Gamma); roii  = find(gamma >= (max(gamma)*TH));
                
%                 T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
%                 J  = zeros(size(LFD, 2), size(Phi_az, 2));
%                 J(roii, :)  = T(roii, :) * Phi_az;
%                 CDR  = currize(J);

%                 rLFD  = LFD(:, roii);
%                 [currycdr] = solve_inverse_problem(Phi_bl, Phi_az, rLFD, [], [], [], ANI, 'rot', 'SBL-SVD', 3);
%                 CDR  = zeros(4, size(LFD, 2)/3, size(Phi_az, 2));
%                 CDR(:, roii(3:3:end)/3, :)  = currycdr;
%                 J  = currycdr(1:3, :, :) .* repmat(currycdr(4, :, :), [3, 1, 1]);
%                 J  = reshape(J, [size(J, 1)*size(J, 2), size(J, 3)]);
%                 figure, imagesc(abs(mypermute(find_cov(J, 'EM'), roii, ANI, 'rot')));           

                rLFD  = LFD(:, roii);
                [~, rGamma] = solve_inverse_problem(Phi_bl, Phi_az, rLFD, [], [], [], ANI, 'rot', 'SBL-SVD', 3);
                Gamma  = sparse(size(LFD, 2), size(LFD, 2)); Gamma(roii, roii)  = rGamma;
                T  = Gamma * LFD.' * (LFD * Gamma * LFD.' + C_bl) ^ -1;
                J = T * Phi_az;
                CDR  = currize(J, [], 'rot');


            %=============================================================%%=============================================================%
            case 'FLCMV-ROI'

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    for j = 0 : 2
                        T(i+j, :) = ((LFD(:, i+j).' * iC_az * LFD(:, i+j)) ^ -1) * LFD(:, i+j).' * iC_az;
                    end
                end

                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ROI'

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                end

                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-DN-ROI'
 
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                iC_bl  = C_bl^-1;
                iC_az  = C_az^-1;
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                    T(i:i+2, :) = T(i:i+2, :) / sqrt(trace((LFD(:, i:i+2).' * iC_bl * LFD(:, i:i+2)) ^ -1));
                end 	 
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-DL-ROI'

                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'EM');
                lambda  = 10 * mean(diag(C_bl));
                iC_az   = (C_az + lambda * eye(size(C_az, 1))) ^-1;
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :) = ((LFD(:, i:i+2).' * iC_az * LFD(:, i:i+2)) ^ -1) * LFD(:, i:i+2).' * iC_az;
                end

                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV-ROI'
                
                gammamat   = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                ExtIDX   = extend_ind(double(IDX), 'rot');
                T  = zeros(size(LFD, 2), size(LFD, 1));
                T(ExtIDX, :)  = mini_rLCMV(C_az, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');  


            %=============================================================%%=============================================================%
            case 'wrLCMV-ROI'
                
                gammamat  = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                ExtIDX   = extend_ind(double(IDX), 'rot');
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, T(ExtIDX, :)]  = mini_rLCMV(C_az, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV-ND-DN-ROI'
                
                gammamat   = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                ExtIDX   = extend_ind(double(IDX), 'rot');
                T     = zeros(size(LFD, 2), size(LFD, 1));
                T_az  = zeros(size(LFD, 2), size(LFD, 1));
                T_bl  = zeros(size(LFD, 2), size(LFD, 1));
                T_az(ExtIDX, :)  = mini_rLCMV(C_az, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                T_bl(ExtIDX, :)  = mini_rLCMV(C_bl, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :) = T_az(i:i+2, :) / sqrt(trace(T_bl(i:i+2, :) * C_bl * T_bl(i:i+2, :).'));
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');                
 

            %=============================================================%%=============================================================%
            case 'wrLCMV-ND-DN-ROI'
                
                gammamat   = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                ExtIDX   = extend_ind(double(IDX), 'rot');
                T     = zeros(size(LFD, 2), size(LFD, 1));
                T_az  = zeros(size(LFD, 2), size(LFD, 1));
                T_bl  = zeros(size(LFD, 2), size(LFD, 1));
                [~, T_az(ExtIDX, :)]  = mini_rLCMV(C_az, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                [~, T_bl(ExtIDX, :)]  = mini_rLCMV(C_bl, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :) = T_az(i:i+2, :) / sqrt(trace(T_bl(i:i+2, :) * C_bl * T_bl(i:i+2, :).'));
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, [], 'rot');  
                
                
            %=============================================================%%=============================================================%
            case 'c-rLCMV-ROI'
                
                gammamat   = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                ExtIDX   = extend_ind(double(IDX), 'rot');
                T      = zeros(size(LFD, 2), size(LFD, 1));
                wT     = zeros(size(LFD, 2), size(LFD, 1));
                T_az   = zeros(size(LFD, 2), size(LFD, 1));
                T_bl   = zeros(size(LFD, 2), size(LFD, 1));
                wT_az  = zeros(size(LFD, 2), size(LFD, 1));
                wT_bl  = zeros(size(LFD, 2), size(LFD, 1));                
                [T_az(ExtIDX, :), wT_az(ExtIDX, :)]  = mini_rLCMV(C_az, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                [T_bl(ExtIDX, :), wT_bl(ExtIDX, :)]  = mini_rLCMV(C_bl, LFD(:, ExtIDX), gammamat(:, :, ExtIDX), 'rot');
                for idx = 1 : length(IDX)
                    i  = 3 * IDX(idx) - 2;
                    T(i:i+2, :)  = T_az(i:i+2, :)  / sqrt(trace(T_bl(i:i+2, :)  * C_bl *  T_bl(i:i+2, :).'));
                    wT(i:i+2, :) = wT_az(i:i+2, :) / sqrt(trace(wT_bl(i:i+2, :) * C_bl * wT_bl(i:i+2, :).'));
                end
                T   = project_onto_pcs( T, C_az, 'Kaiser', 1);
                wT  = project_onto_pcs(wT, C_az, 'Kaiser', 1);
                J  = T_az  * Phi_az; CDR(:, :, :, 1)  = currize(J, [], 'rot');
                J  = wT_az * Phi_az; CDR(:, :, :, 2)  = currize(J, [], 'rot');
                J  = T  * Phi_az   ; CDR(:, :, :, 3)  = currize(J, [], 'rot');
                J  = wT * Phi_az   ; CDR(:, :, :, 4)  = currize(J, [], 'rot'); 
                
                
            %=============================================================%%=============================================================%
            case 'c-rLCMV'
                
                gammamat   = extend_gammamat(varargin{4}, size(LFD, 2)/size(varargin{4}, 3));
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                T      = zeros(size(LFD, 2), size(LFD, 1));
                wT     = zeros(size(LFD, 2), size(LFD, 1));
                T_az   = zeros(size(LFD, 2), size(LFD, 1));
                T_bl   = zeros(size(LFD, 2), size(LFD, 1));
                wT_az  = zeros(size(LFD, 2), size(LFD, 1));
                wT_bl  = zeros(size(LFD, 2), size(LFD, 1));                
                [T_az, wT_az]  = mini_rLCMV(C_az, LFD, gammamat, 'rot');
                [T_bl, wT_bl]  = mini_rLCMV(C_bl, LFD, gammamat, 'rot');
                for i = 1 : 3 : size(LFD, 2)
                    T(i:i+2, :)  = T_az(i:i+2, :)  / sqrt(trace(T_bl(i:i+2, :)  * C_bl *  T_bl(i:i+2, :).'));
                    wT(i:i+2, :) = wT_az(i:i+2, :) / sqrt(trace(wT_bl(i:i+2, :) * C_bl * wT_bl(i:i+2, :).'));
                end
                T   = project_onto_pcs( T, C_az, 'Kaiser', 1);
                wT  = project_onto_pcs(wT, C_az, 'Kaiser', 1);
                J  = T_az  * Phi_az; CDR(:, :, :, 1)  = currize(J, [], 'rot');
                J  = wT_az * Phi_az; CDR(:, :, :, 2)  = currize(J, [], 'rot');
                J  = T  * Phi_az   ; CDR(:, :, :, 3)  = currize(J, [], 'rot');
                J  = wT * Phi_az   ; CDR(:, :, :, 4)  = currize(J, [], 'rot'); 
                
                          
             %=============================================================%%=============================================================%
                        

        end    
            
        
%=========================================================================%%=========================================================================%
    case 'fxd'

        switch cdr_type

            %=============================================================%%=============================================================%
            case 'LS'
                
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, S, V]  = svd(A, 'econ');
                J  = V * (S^-1) * U.' * b;
                CDR  = currize(J, ORI, 'fxd');
                
                
            %=============================================================%%=============================================================%
            case 'LS-ROI'
                                
                C_bl = find_cov(Phi_bl, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', zeros(1, nnd));
                LFD  = LFD(:, IDX);
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, S, V]  = svd(A, 'econ');
                J  = zeros(size(ORI, 2), size(b, 2));
                J(IDX, :)  = V * (S^-1) * U.' * b;
                CDR  = currize(J, ORI, 'fxd');
            

            %=============================================================%%=============================================================%
            case 'WMN-DCP'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                w  = norms(A) .^ 2;
                A  = A .* repmat((w.^(-.5)), [size(A, 1), 1]);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta);  
                lambda  = lambda ^ 2;
                T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
                J  = repmat((w.^(-.5)).', [1, size(b, 2)]) .* (T * b);
                CDR  = currize(J, ORI, 'fxd');
                
                                                                    
            %=============================================================%%=============================================================%
            case 'sLORETA-DCP'
 
                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta); 
                lambda  = lambda ^ 2;
                F1  = repmat(s ./ (s.^2 + lambda), [1, size(b, 2)]);
                F2  = repmat(s.^2 ./ (s.^2 + lambda), [1, size(b, 2)]);
                J  = V * (F1 .* (U.' * b));
                C  = V.^2 * F2;
                CDR  = currize((J.^2) ./ C, ORI, 'fxd');

                
            %=============================================================%%=============================================================%
            case 'sLORETA-DCP-ROI'
 
                C_bl  = find_cov(Phi_bl, 'SC');                
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                [A, b]  = mywhitening(LFD(:, IDX), Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                [~, azp_ind]  = max(var(b));
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                [~, lambda] = discrep(U, s, V, b(:, azp_ind), delta); 
                lambda  = lambda ^ 2;
                F1  = repmat(s ./ (s.^2 + lambda), [1, size(b, 2)]);
                F2  = repmat(s.^2 ./ (s.^2 + lambda), [1, size(b, 2)]);
                J  = V * (F1 .* (U.' * b));
                C  = V.^2 * F2;
                CDR  = zeros(4, size(LFD, 2), size(b, 2));
                CDR(:, IDX, :)  = currize((J.^2) ./ C, ORI(:, IDX), 'fxd');
                
                
            %=============================================================%%=============================================================%    
            case 'sLORETA-DCP-TC'

                C_bl  = find_cov(Phi_bl, 'SC');
                [A, b]  = mywhitening(LFD, Phi_az, C_bl);
                [U, s, V]  = csvd(A);
                delta  = sqrt(chi2inv(0.975, size(b, 1)));
                F1  = zeros(size(Phi_az));
                F2  = zeros(size(Phi_az));
                for i = 1 : size(Phi_az, 2)
                    [~, lambda] = discrep(U, s, V, b(:, i), delta); 
                    lambda  = lambda ^ 2;
                    F1(:, i)  = s ./ (s.^2 + lambda);
                    F2(:, i)  = s.^2 ./ (s.^2 + lambda);
                end
                J  = V * (F1 .* (U.' * b));
                C  = V.^2 * F2;
                CDR  = currize((J.^2) ./ C, ORI, 'fxd');
                
                
            %=============================================================%%=============================================================%
            case 'LCMV'

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                P  = LFD.' * iC_az;
                T  = repmat(dot(P, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P;
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');                
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND'

                iC_bl  = find_cov(Phi_bl, 'SC')^-1;
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                P_az  = LFD.' * iC_az;
                P_bl  = LFD.' * iC_bl;
                T  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T  = T ./ repmat(dot(P_bl, LFD.', 2).^-0.5, [1, size(iC_az, 2)]);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 

                
            %=============================================================%%=============================================================%
            case 'LCMV-DN'
 
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                P_az  = LFD.' * iC_az;
                T  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-DN'
 
                C_bl  = find_cov(Phi_bl, 'SC'); iC_bl  = C_bl^-1;
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                P_az  = LFD.' * iC_az;
                P_bl  = LFD.' * iC_bl;
                T  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T  = T ./ repmat(dot(P_bl, LFD.', 2).^-0.5, [1, size(iC_az, 2)]);
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');  
                
                
            %=============================================================%%=============================================================%
            case 'wLCMV'

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                LFD  = LFD ./ repmat(norms(LFD), [size(LFD, 1), 1]);
                P  = LFD.' * iC_az;
                T  = repmat(dot(P, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P;
   
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'wLCMV-DN'
 
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                LFD  = LFD ./ repmat(norms(LFD), [size(LFD, 1), 1]);
                P_az  = LFD.' * iC_az;
                T  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ROI'
                    
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                LFD  = LFD(:, IDX);
                P  = LFD.' * iC_az;
                T(IDX, :) = repmat(dot(P, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P;
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 

                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-ROI'

                iC_bl  = find_cov(Phi_bl, 'SC')^-1;
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                LFD  = LFD(:, IDX);
                P_az  = LFD.' * iC_az;
                P_bl  = LFD.' * iC_bl;
                T(IDX, :)  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T(IDX, :)  = T(IDX, :) ./ repmat(dot(P_bl, LFD.', 2).^-0.5, [1, size(iC_az, 2)]);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                

            %=============================================================%%=============================================================%
            case 'LCMV-DN-ROI'
 
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                LFD  = LFD(:, IDX);
                P_az  = LFD.' * iC_az;
                T(IDX, :)  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-DN-ROI'
 
                C_bl  = find_cov(Phi_bl, 'SC'); iC_bl  = C_bl^-1;
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                LFD  = LFD(:, IDX);
                P_az  = LFD.' * iC_az;
                P_bl  = LFD.' * iC_bl;
                T(IDX, :)  = repmat(dot(P_az, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P_az;
                T(IDX, :)  = T(IDX, :) ./ repmat(dot(P_bl, LFD.', 2).^-0.5, [1, size(iC_az, 2)]);
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                    
            %=============================================================%%=============================================================%
            case 'wLCMV-ROI'
                    
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                LFD  = LFD ./ repmat(norms(LFD), [size(LFD, 1), 1]);
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                LFD  = LFD(:, IDX);
                P  = LFD.' * iC_az;
                T(IDX, :) = repmat(dot(P, LFD.', 2).^-1, [1, size(iC_az, 2)]) .* P;
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');               
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV'
                
                gammamat   = varargin{4}; % leadfield noise covariance
                C_az  = find_cov(Phi_az, 'SC');
                T  = mini_rLCMV(C_az, LFD, gammamat, 'fxd');
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV-ROI'
                
                gammamat   = varargin{4}; % leadfield noise covariance
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                T(IDX, :) = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 

                
            %=============================================================%%=============================================================%
            case 'wrLCMV-ROI'
                
                gammamat   = varargin{4}; % leadfield noise covariance
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, T(IDX, :)]  = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');
                

            %=============================================================%%=============================================================%
            case 'rLCMV-ND-ROI'

                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T     = zeros(size(LFD, 2), size(LFD, 1));
                T_az  = zeros(size(LFD, 2), size(LFD, 1));
                T_bl  = zeros(size(LFD, 2), size(LFD, 1));
                T_az(IDX, :)  = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                T_bl(IDX, :)  = mini_rLCMV(C_bl, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                for i_idx = 1 : length(IDX)
                    T(IDX(i_idx), :)  = T_az(IDX(i_idx), :) / sqrt(T_bl(IDX(i_idx), :) * C_bl * T_bl(IDX(i_idx), :).');
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');  
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV-DN-ROI'
                
                gammamat   = varargin{4}; % leadfield noise covariance
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                T(IDX, :) = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');
                
                
            %=============================================================%%=============================================================%
            case 'rLCMV-ND-DN-ROI'

                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T     = zeros(size(LFD, 2), size(LFD, 1));
                T_az  = zeros(size(LFD, 2), size(LFD, 1));
                T_bl  = zeros(size(LFD, 2), size(LFD, 1));
                T_az(IDX, :)  = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                T_bl(IDX, :)  = mini_rLCMV(C_bl, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                for i_idx = 1 : length(IDX)
                    T(IDX(i_idx), :)  = T_az(IDX(i_idx), :) / sqrt(T_bl(IDX(i_idx), :) * C_bl * T_bl(IDX(i_idx), :).');
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');  
                
                
            %=============================================================%%=============================================================%
            case 'wrLCMV-ND-DN-ROI'

                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T     = zeros(size(LFD, 2), size(LFD, 1));
                T_az  = zeros(size(LFD, 2), size(LFD, 1));
                T_bl  = zeros(size(LFD, 2), size(LFD, 1));
                [~, T_az(IDX, :)]  = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                [~, T_bl(IDX, :)]  = mini_rLCMV(C_bl, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                for i_idx = 1 : length(IDX)
                    T(IDX(i_idx), :)  = T_az(IDX(i_idx), :) / sqrt(T_bl(IDX(i_idx), :) * C_bl * T_bl(IDX(i_idx), :).');
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd');
                
                
                
            %=============================================================%%=============================================================%
            case 'c-rLCMV-ROI'
            
                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T      = zeros(size(LFD, 2), size(LFD, 1));
                wT     = zeros(size(LFD, 2), size(LFD, 1));
                T_az   = zeros(size(LFD, 2), size(LFD, 1));
                T_bl   = zeros(size(LFD, 2), size(LFD, 1));
                wT_az  = zeros(size(LFD, 2), size(LFD, 1));
                wT_bl  = zeros(size(LFD, 2), size(LFD, 1));                
                [T_az(IDX, :), wT_az(IDX, :)]  = mini_rLCMV(C_az, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');
                [T_bl(IDX, :), wT_bl(IDX, :)]  = mini_rLCMV(C_bl, LFD(:, IDX), gammamat(:, :, IDX), 'fxd');              
                for i_idx = 1 : length(IDX)
                    T(IDX(i_idx), :)   = T_az(IDX(i_idx), :)  / sqrt(T_bl(IDX(i_idx), :)  * C_bl *  T_bl(IDX(i_idx), :).');
                    wT(IDX(i_idx), :)  = wT_az(IDX(i_idx), :) / sqrt(wT_bl(IDX(i_idx), :) * C_bl * wT_bl(IDX(i_idx), :).');
                end
                T   = project_onto_pcs( T, C_az, 'Kaiser', 1);
                wT  = project_onto_pcs(wT, C_az, 'Kaiser', 1);
                J  = T_az  * Phi_az; CDR(:, :, :, 1)  = currize(J, ORI, 'fxd');
                J  = wT_az * Phi_az; CDR(:, :, :, 2)  = currize(J, ORI, 'fxd'); 
                J  = T  * Phi_az   ; CDR(:, :, :, 3)  = currize(J, ORI, 'fxd');
                J  = wT * Phi_az   ; CDR(:, :, :, 4)  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'c-rLCMV'
            
                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');
                T      = zeros(size(LFD, 2), size(LFD, 1));
                wT     = zeros(size(LFD, 2), size(LFD, 1));
                T_az   = zeros(size(LFD, 2), size(LFD, 1));
                T_bl   = zeros(size(LFD, 2), size(LFD, 1));
                wT_az  = zeros(size(LFD, 2), size(LFD, 1));
                wT_bl  = zeros(size(LFD, 2), size(LFD, 1));                
                [T_az, wT_az]  = mini_rLCMV_fxd(C_az, LFD, gammamat);
                [T_bl, wT_bl]  = mini_rLCMV_fxd(C_bl, LFD, gammamat);              
                T   = T_az   ./ repmat(sqrt(dot(T_bl,  T_bl*C_bl, 2)),  [1, size(LFD, 1)]);
                wT  = wT_az  ./ repmat(sqrt(dot(wT_bl, wT_bl*C_bl, 2)), [1, size(LFD, 1)]);
                T   = project_onto_pcs( T, C_az, 'Kaiser', 1);
                wT  = project_onto_pcs(wT, C_az, 'Kaiser', 1);
                J  = T_az  * Phi_az; CDR(:, :, :, 1)  = currize(J, ORI, 'fxd');
                J  = wT_az * Phi_az; CDR(:, :, :, 2)  = currize(J, ORI, 'fxd'); 
                J  = T  * Phi_az   ; CDR(:, :, :, 3)  = currize(J, ORI, 'fxd');
                J  = wT * Phi_az   ; CDR(:, :, :, 4)  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-SBP' % Shahbazpanahi

                gammamat  = varargin{4};
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : size(LFD, 2)
                    sigma    = mean(diag(gammamat(:, :, i)));
                    epsilon  = sqrt(size(LFD, 1)) * (sigma^2);
                    T(i, :)  = mini_sbp(LFD(:, i), iC_az, epsilon);
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-SBP' % Shahbazpanahi

                gammamat  = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC'); iC_bl  = C_bl^-1;
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                W  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : size(LFD, 2)
                    sigma    = mean(diag(gammamat(:, :, i)));
                    epsilon  = sqrt(size(LFD, 1)) * (sigma^2);
                    T(i, :)  = mini_sbp(LFD(:, i), iC_az, epsilon);
                    W(i, :)  = mini_sbp(LFD(:, i), iC_bl, epsilon);
                    T(i, :)  = T(i, :) / sqrt(W(i, :) * C_bl * W(i, :).');
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-ND-DN-SBP' % Shahbazpanahi

                gammamat  = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC'); iC_bl  = C_bl^-1;
                C_az  = find_cov(Phi_az, 'SC'); iC_az  = C_az^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                W  = zeros(size(LFD, 2), size(LFD, 1));
                for i = 1 : size(LFD, 2)
                    sigma   = mean(diag(gammamat(:, :, i)));
                    epsilon  = sqrt(size(LFD, 1)) * (sigma^2);
                    T(i, :)  = mini_sbp(LFD(:, i), iC_az, epsilon);
                    W(i, :)  = mini_sbp(LFD(:, i), iC_bl, epsilon);
                    T(i, :)  = T(i, :) / sqrt(W(i, :) * C_bl * W(i, :).');
                end
                T  = project_onto_pcs(T, C_az, 'Kaiser', 1);
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-SBP-ROI' % Shahbazpanahi

                gammamat  = varargin{4};
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                for i_idx = 1 : length(IDX)
                    sigma   = mean(diag(gammamat(:, :, IDX(i_idx))));
                    epsilon  = sqrt(size(LFD, 1)) * (sigma^2);
                    T(IDX(i_idx), :)  = mini_sbp(LFD(:, IDX(i_idx)), iC_az, epsilon);
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'LCMV-SBP-Test' % Shahbazpanahi

                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                
                for i_idx = 1 : length(IDX)
                    [LFD(:, IDX(i_idx)), epsilon]  = my_temporary_add_noise(LFD(:, IDX(i_idx)), 20);
                    T(IDX(i_idx), :)  = mini_sbp(LFD(:, IDX(i_idx)), iC_az, epsilon);
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
              
                
            %=============================================================%%=============================================================%
            case 'LCMV-SBP-ROI-EST' % Shahbazpanahi

                gammamat  = varargin{4};
                cleanlfd  = varargin{4};
                iC_az  = find_cov(Phi_az, 'SC')^-1;
                T  = zeros(size(LFD, 2), size(LFD, 1));
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                for i_idx = 1 : length(IDX)
                    sigma    = mean(diag(gammamat(:, :, IDX(i_idx))));
                    epsilon  = estimate_epsilon(cleanlfd(:, IDX(i_idx)), sigma, 10000);
                    T(IDX(i_idx), :)  = mini_sbp(LFD(:, IDX(i_idx)), iC_az, epsilon);
                end
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'myrLCMV-ROI'

                gammamat   = varargin{4}; % leadfield noise covariance
                C_az  = find_cov(Phi_az, 'SC');
                [~, IDX] = find_vicinity(varargin{1}.sdl, varargin{2}, TRI, 'area', ROIA*ones(1, nnd));
                T  = zeros(size(LFD, 2), size(LFD, 1));
                T(IDX, :) = mini_myrLCMV(C_az, LFD(:, IDX), gammamat(IDX), 'fxd');
                J  = T * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
         
            %=============================================================%%=============================================================%
            case 'RMVB'
            
                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');              
                [T_az, wT_az]  = mini_rLCMV_fxd(C_az, LFD, gammamat);
                [T_bl, wT_bl]  = mini_rLCMV_fxd(C_bl, LFD, gammamat);   
                J  = wT_az * Phi_az; 
                CDR  = currize(J, ORI, 'fxd'); 
                
                
            %=============================================================%%=============================================================%
            case 'RMVB-ND-DN'

                gammamat   = varargin{4};
                C_bl  = find_cov(Phi_bl, 'SC');
                C_az  = find_cov(Phi_az, 'SC');               
                [T_az, wT_az]  = mini_rLCMV_fxd(C_az, LFD, gammamat);
                [T_bl, wT_bl]  = mini_rLCMV_fxd(C_bl, LFD, gammamat);              
                T   = T_az   ./ repmat(sqrt(dot(T_bl,  T_bl*C_bl, 2)),  [1, size(LFD, 1)]);
                wT  = wT_az  ./ repmat(sqrt(dot(wT_bl, wT_bl*C_bl, 2)), [1, size(LFD, 1)]);
                T   = project_onto_pcs( T, C_az, 'Kaiser', 1);
                wT  = project_onto_pcs(wT, C_az, 'Kaiser', 1);
                J  = wT * Phi_az;
                CDR  = currize(J, ORI, 'fxd'); 
                
                
        end
 
        
%=========================================================================%%=========================================================================%
    otherwise
        display('*************** Weird Request ***************')
        
        
end


end



%%

function epsilon = estimate_epsilon(lfd, sigma, N)
rng(1)
M = size(lfd, 1);
err = zeros(N, 1);
for i = 1 : N
    tmp = lfd + sigma * randn(M, 1);
    err(i)  = norm(tmp*tmp.' - lfd*lfd.', 'fro');
%     err(:, :, i)  = tmp*tmp.' - lfd*lfd.';
end
% epsilon  = mean(err) + sqrt(chi2inv(0.975, 1)) * std(err, 1);
epsilon  = mean(err);
figure, histogram(err), hold on, stem(epsilon, 400)
end


function T = mini_sbp(lfd, iC, epsilon)
if (epsilon >= max(eig(lfd * lfd.')))
    display('warning!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
end
Rs  = lfd * lfd.' - epsilon * eye(size(lfd, 1));
[a, ~, ~]  = svds(iC * Rs, 1);
% T  = a;
T  = ((a.' * Rs * a)^-0.5) * a;
T  = T / sign(T.' * lfd);
end


function [T, wT] = mini_rLCMV_fxd(C, LFD, gammamat)

[m, n]  = size(LFD); T  = zeros(n, m); wT  = zeros(n, m);

%%
parfor i = 1 : n
    A  = gammamat(:, :, i);  lfd  = LFD(:, i);
    w  = my_little_cvx(C, lfd, A);
%     w  = my_little_opt(C, lfd, A); % faster version
    T(i, :)  = w.';
    wT(i, :) = ((T(i, :)*lfd)^-1) * T(i, :);
end

end

function w = my_little_opt(R, c, A)
% For the theory see R. G. Lorenz and S. P. Boyd, “Robust minimum variance beamforming,” IEEE Trans. Signal Process., vol. 53, no. 5, pp. 1684–1696, May 2005. 
% (https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1420809) for the

    Q         = A*A.'- c*c.';    
    [~, iRh]  = myx2factrzn(R);
    Qt        = iRh * Q * iRh.';
    [V, G]    = eig(Qt); 
    g         = diag(G);
    ct        = V.' * iRh * c;
    fun       = @(lambda)...
                (lambda^2)*sum((ct.^2).*g./((1+lambda*g).^2))-(2*lambda)*sum((ct.^2)./(1+lambda*g))-1;
    cj        = ct(g<0); 
    gj        = g(g<0); 
    lambda0   = -(1+abs(cj)/sqrt(gj+cj^2))/gj;
    lambda    = fzero(fun, lambda0);
    w         = -lambda * ((R + lambda * Q)^-1) * c; 
    
%     t=lambda0:0.1*lambda0:2000*lambda0;
%     for tt=1:length(t), ft(tt)=fun(t(tt)); end
%     figure, plot(t, ft), hold on, stem(lambda, 1, 'b')
end


% function w = my_little_opt(R, c, A)
%     Q         = A*A.'- c*c.';    
%     [~, iRh]  = myx2factrzn(R);
%     Qt        = iRh * Q * iRh.';
%     [V, G]    = eig(Qt); 
%     g         = diag(G);
%     ct        = V.' * iRh * c;
% %     fun       = @(lambda)...
% %         [(lambda^2)*sum((ct.^2).*g./((1+lambda*g).^2))-(2*lambda)*sum((ct.^2)./(1+lambda*g))-1, ...
% %         -2*sum((ct.^2)./((1+lambda*g).^3))];
%     fun1      = @(lambda)...
%         (lambda^2)*sum((ct.^2).*g./((1+lambda*g).^2))-(2*lambda)*sum((ct.^2)./(1+lambda*g))-1;
%     cj        = ct(g<0); 
%     gj        = g(g<0); 
%     lambda0   = -(1+abs(cj)/sqrt(gj+cj^2))/gj;
% %     [lambda1, resnorm, ~, exitflag, ~, ~] = newtonraphson(fun, lambda0);
%     [lambda, fval, exitflag, ~] = fzero(fun1, lambda0);
% %     fval
% %     exitflag
%     w         = -lambda * ((R + lambda * Q)^-1) * c; 
%     
% %     t=lambda0:0.1*lambda0:2*lambda0;
% %     for tt=1:length(t), ft(tt)=fun1(t(tt)); end
% %     figure, plot(t, ft), hold on, stem(lambda, 1, 'b'), hold on, stem(lambda1, 1, 'r')
% end


function w = my_little_cvx(C, lfd, A)
% lfd.'*((A*A.')^-1)*lfd
iter = 1;
while(1)
    cvx_begin quiet
        variable w(size(C, 1), 1)
        minimize( w.' * C * w )
        subject to
            w.' * lfd - norms(A.' * w) >= 1;
    cvx_end
    if(~isinf(w) & ~isnan(w)), break, end
    iter = iter + 1;
    A = A * 0.75;
end
end

function [T, wT]  = mini_rLCMV(C, LFD, gammamat, lfd_type)

if(strcmp(lfd_type, 'rot')), d = 3; elseif(strcmp(lfd_type, 'fxd')), d = 1; end
m  = size(LFD, 1);
C  = mysparseblockdiag(repmat(C, [d, 1]));
T  = zeros(size(LFD, 2), size(LFD, 1)); wT = zeros(size(LFD, 2), size(LFD, 1));

%%
% start_point  = tic;
% h = waitbar(0, 'Preparing to solve RMVB problem...');
for i = 1 : d : size(LFD, 2)

    if(~gammamat(:, :, i:i+d-1)), error('Warning: point uncertainty region!!!'), end
    A  = gammamat(:, :, i:i+d-1);  lfd  = LFD(:, i:i+d-1);
    delta  = 0.025;  iter  = 0; iter_max  = 10;

    while(iter < 10)
        cvx_begin quiet
            variable W(size(C, 1), 1)
            minimize( W.' * C * W )
            subject to
                for j = 1 : d
                    wj   = W((j-1)*m+1:j*m);
                    for k = 1 : d
                        if(k == j)
                            wj.' * lfd(:, j) - norms(A(:, :, j).' * wj) >= 1;
                        else
                            abs(wj.' * lfd(:, k)) + norms(A(:, :, k).' * wj) <= delta;
                        end
                    end                        
                end
        cvx_end
        if(~sum(isnan(W))), break; end
        display('*************warning: iterative procedure!!!*************')
        delta  = delta + 0.025;
        iter = iter + 1
    end

    T(i:i+d-1, :)  = reshape(W, [], d).';
    wT(i:i+d-1, :) = ((T(i:i+d-1, :)*lfd)^-1) * T(i:i+d-1, :);
%     mini_rLCMV_display(lfd, AT.', W) 
%     percentage  = i/size(LFD, 2);
%     waitbar(percentage, h, ['Solving wrLCMV problem (' datestr((1-percentage)/percentage*toc(start_point)/60/60/24, 'HH:MM:SS') ' left...']); 
end
% close(h)

end

function mini_rLCMV_display(lfd, A, W)
    d = size(lfd, 2);  
    W  = reshape(W, [], d);
    mat_worst_case = zeros(d);
    mat_bestt_case = zeros(d);
    for h = 1 : d
        for p = 1 : d
            U = squeeze(A(:, :, p)).' * W(:, h); U = U/norms(U);
            E(:, p) = -squeeze(A(:, :, p)) * U; 
        end
        mat_worst_case(h, :) = W(:, h).' * (lfd + E);
        mat_bestt_case(h, :) = W(:, h).' * (lfd - E);
    end
    display('============================================================')
    mat_worst_case
    mat_bestt_case
    mat_centr_case  = (W.' * lfd)
    eig_value_case  = eig(mat_centr_case - eye(d));
    display('============================================================')
end


function T = project_onto_pcs(T, C, options, param)
[U, S, ~] = csvd(C); 
switch options
    case 'user'
        est_nnd  = param; 
        
    case 'Kaiser'        
        est_nnd  = find(S/mean(S) >= param);  
     
    case 'screetest'    
        
    case 'minvariance'
        est_nnd  = find(cumsum(S)/sum(S) <= param);
        
end        
if(isempty(est_nnd)), est_nnd = 1; else, est_nnd  = est_nnd(end); end
T  = U(:, 1:est_nnd) * U(:, 1:est_nnd).' * T.';  
T  = T.';
end



function gammamat   = extend_gammamat(gammamat, d)

m  = size(gammamat, 1);
gammamat  = repmat(gammamat, [1, d, 1]);
gammamat  = reshape(gammamat, m, m, []);

end


function epi  = find_nvoxel_interface(epl, curryloc)

epi = cell(1, length(epl));
for i = 1 : length(epl)
    [~, epi{i}]  = find_nvoxel(epl{i}, curryloc);
end

end



function [A, b, iCh] = mywhitening(A, b, C)

[U, S, ~]  = svd(C^-1);
iCh  = sqrt(S) * U.';
A  =  iCh * A;
b  =  iCh * b;

end

function CDR  = currize(J, ori, type, varargin)

switch type
    case 'rot'
        J  = reshape(J, [3, size(J, 1)/3, size(J,2)]);
        CDR  = zeros(4, size(J, 2), size(J, 3));
        CDR(1:3, :, :)  = J ./ repmat(norms(J, 2, 1), [3, 1, 1]);
        CDR(4, :, :)  = norms(J, 2, 1);
        CDR(isnan(CDR))  = 0;
        
    case 'fxd'
        CDR  = zeros(4, size(J, 1), size(J, 2));
        CDR(1:3, :, :) = repmat(ori, [1, 1, size(J, 2)]);
        CDR(4, :, :)  = J;
end

if(~isempty(varargin))
    CDR  = repmat(CDR, [1, 1, varargin{:}]);
end
end


function [X, iX]  = myx2factrzn(C)

[U, S, ~]  = svd(C);
X    = U * (S .^ 0.5) * U.';
iX   = U * diag(diag(S) .^ -0.5) * U.';

end


function X  = myxxtfactrzn(C, varargin)

[U, S, ~]  = svd(C);
if(isempty(varargin))
    r  = nnz(diag(S));
else
    r  = varargin{:};
end
S  = sqrt(S);
X  = U(:, 1:r) * S(1:r, 1:r);

end



function Gamma  = mysparseblockdiag(gamma)
% input should be vertical
[m, n]  = size(gamma);
I  = repmat((1:m).', [1, n]);
J  = reshape(1:m, [n, m/n]);
J  = repmat(J, [n, 1]);
J  = reshape(J, [n,  m]);
J  = J.';
Gamma  = sparse(I, J, gamma, m, m);

end

function gamma  = myblockexpansion(gamma, d)

gamma  = reshape(gamma, [1, length(gamma)]);
gamma  = repmat(gamma, [d, 1]);
gamma  = gamma(:); 

end

function gamma  = myblockaverage(gamma, d)

gsize  = size(gamma);
gamma  = reshape(gamma, [1, length(gamma)]);
gamma  = reshape(gamma, d, []);
gamma  = mean(gamma, 1);
gamma  = repmat(gamma, [d, 1]);
gamma  = reshape(gamma, gsize); 

end

function output  = mygblockaverage(gamma, ANI)

output = zeros(size(gamma));
for i = 1 : length(ANI)
    ind  = extend_ind(double(ANI{i}), 'rot');
    output(ind)  = mean(gamma(ind));
%     output(ind)  = max(gamma(ind), [], 1);
end

end

function gamma  = myggblockaverage(gamma, TRI)

for i = 1 : size(TRI, 2)
    ind  = extend_ind(double(TRI(:, i)).', 'rot');
    gamma(ind([1:3 4:6]))  = mean(gamma(ind([1:3 4:6])));
    gamma(ind([4:6 7:9]))  = mean(gamma(ind([4:6 7:9])));
    gamma(ind([1:3 7:9]))  = mean(gamma(ind([1:3 7:9])));
end

end

function [J, lambda] = mywmn(A, b, w)

A  = A * mysparseblockdiag(w.^0.5);
[U, s, V]  = csvd(A);
[~, azp_ind]  = max(var(b));
delta  = sqrt(chi2inv(0.975, size(b, 1)));
[~, lambda] = discrep(U, s, V, b(:, azp_ind), delta);  
lambda  = lambda ^ 2;
T  = A.' * (A * A.' + lambda * eye(size(b, 1))) ^ -1;
J  = mysparseblockdiag(w.^0.5) * (T * b);   

end

function C  = mypermute(C, I, ANI, lfdtype)
idx  = [];
for i = 1 : length(ANI)
    ind  = extend_ind(double(ANI{i}), lfdtype);
    idx  = [idx; find(ismember(I, ind))];
end
C  = C(idx, :);
C  = C(:, idx);
end

function gamma  = solve_sparse_problem(gamma, x, a, b, NHD, d)

if(numel(b) == 1), b = b * ones(size(a)); end

for i = 1 : length(gamma)
    P  = sort(gamma(NHD{i}), 'ascend'); G  = [0; P; inf]; % potential values
    for j = 1 : length(G) - 1 
        g  = x(i) / sqrt(a(i) + b(i)*sum(sign((G(j)+G(j+1))/2 - gamma(NHD{i}))));
        if(isreal(g) && g>G(j) && g<G(j+1)), P  = [P; g]; end
    end
    F  = (x(i)^2)./P + a(i)*P + b(i)*sum(abs(repmat(P, [1, length(NHD{i})])-repmat(gamma(NHD{i}).', [length(P), 1])), 2);
    [~, ind] = min(F);
    if(isempty(ind))
        i
        NHD{i}
        gamma(NHD{i})
        P
        F
        ind
        P(ind)
    end
    gamma(i) = P(ind);
end

gamma  = myblockexpansion(gamma, d);

end

function [fval, val1, val2]  = find_fval(NHD, alpha, lambda, option, varargin)

switch option 
    case 'REG' % regular
        gamma = varargin{1};
        C_az  = varargin{2};
        S_az  = varargin{3};
        iS_az = varargin{4};
        val1  = trace(C_az * iS_az) + sum(log(eig(S_az)));

    case 'CVX' % convexified cost function
        gamma = varargin{1};
        X  = varargin{2};
        Z  = varargin{3}; 
        d  = varargin{4};
        val1  = 0;
        for i = 1 : d : length(gamma)
            val1  = val1 + (norm(X(i:i+2, :), 'fro')^2)/gamma(i) + sum(Z(i:i+2, :).'*gamma(i));
        end
        gamma = gamma(1:d:end);
end

variation_cost = 0;
for i = 1 : length(NHD)
    variation_cost = variation_cost + sum(abs(gamma(i)-gamma(NHD{i})));
end
val2  = alpha*sqrt(3)*sum(gamma)+(1-alpha)*sqrt(3)*variation_cost/2;
fval  = val1 + lambda * val2;

end

function [ws, wv]  = update_weights(gamma, NHD, d)

ws  = abs(gamma);
wv  = zeros(size(gamma));
for i = 1 : length(NHD)
    wv(i) = sum(abs(gamma(i)-gamma(NHD{i}))) / 2;
end

epss  = 0.1 *  max(ws);
epsv  = 0.1 *  max(wv);

ws  = myblockexpansion(epss ./ (ws + epss), d);
wv  = myblockexpansion(epsv ./ (wv + epsv), d);

end




function iC = my3by3inv(C)

% This function efficiently computes the inverse of a series of 3 by 3
% matrices in the form of either C_{3x3x:} or C_{:, 3}
% The formulas are based on https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_.C3.97_3_matrices

%%

if(ndims(C) == 3)

    x1  = C(:, 1, :);
    x2  = C(:, 2, :);
    x3  = C(:, 3, :);
    
    x2Xx3  = permute(cross(x2, x3, 1), [2, 1, 3]);
    x3Xx1  = permute(cross(x3, x1, 1), [2, 1, 3]);
    x1Xx2  = permute(cross(x1, x2, 1), [2, 1, 3]);
    
    DET  = dot(x1, cross(x2, x3, 1), 1);
    DET  = repmat(DET, [3, 3, 1]);
    
    iC  = cat(1, x2Xx3, x3Xx1, x1Xx2) ./ DET;
    
else
    
    C  = reshape(C.', [3, 3, size(C, 1)/3]);
    C  = permute(C, [2, 1, 3]);
    
    x1  = C(:, 1, :);
    x2  = C(:, 2, :);
    x3  = C(:, 3, :);
    
    x2Xx3  = permute(cross(x2, x3, 1), [2, 1, 3]);
    x3Xx1  = permute(cross(x3, x1, 1), [2, 1, 3]);
    x1Xx2  = permute(cross(x1, x2, 1), [2, 1, 3]);
    
    DET  = dot(x1, cross(x2, x3, 1), 1);
    DET  = repmat(DET, [3, 3, 1]);
    
    iC  = cat(1, x2Xx3, x3Xx1, x1Xx2) ./ DET; 
    iC  = permute(iC, [2, 1, 3]);
    iC  = reshape(iC, [3, 3*size(iC, 3)]).';
    
end

end


function [B, L]  = revise_grp(A, ANI, lfd_type, revision_type)

if(strcmp(lfd_type, 'rot')), d = 3; elseif(strcmp(lfd_type, 'fxd')), d = 1; end

switch revision_type
    case 'Forward'
        B  = [];
        L  = zeros(1, length(ANI));
        for i = 1 : length(ANI)
            B  = [B A(:, extend_ind(double(ANI{i}), lfd_type))];
            L(i)  = d * length(ANI{i}); 
        end
    case 'Inverse'
        L  = 0;
        B  = zeros(size(A));
        for i = 1 : length(ANI) 
            B(extend_ind(double(ANI{i}), lfd_type), :) = A(extend_ind(L+1:L+length(ANI{i}), lfd_type), :);
            L  = L + length(ANI{i});
        end
end

end


function J = glasso_cvx(A, b, delta, ANI, costfun_type, alpha, lfd_type)

if(strcmp(lfd_type, 'rot')), d = 3; 
elseif(strcmp(lfd_type, 'fxd')), d = 1; end

switch costfun_type
    case 'regular'
        cvx_begin
            variable J(size(A, 2), size(b, 2))
            z = 0;
            for ind = 1 : length(ANI)
                z = z + norm(J(extend_ind(double(ANI{ind}), lfd_type), :), 'fro') / sqrt(d * length(ANI{ind}));
            end
            minimize( z )
            subject to
                norm(A * J - b, 'fro') <= sqrt(size(b, 2)) * delta;
        cvx_end
        
    case 'sparse GL'
        cvx_begin
            variable J(size(A, 2), size(b, 2))
            z1 = 0;
%             for ind = 1 : length(ANI)
%                 z1 = z1 + (norm(J(extend_ind(double(ANI{ind}), lfd_type), :), 'fro') / sqrt(d * length(ANI{ind})));
%             end
            z2 = 0;
            for ind = 1 : size(A, 2)/d
                z1 = z1 + norm(J(d*(ind-1)+1:d*ind, :), 2) / sqrt(3);
            end
            z  = alpha * z1 + (1-alpha)*z2;
            minimize( z )
            subject to
                norm(A * J - b, 'fro') <= sqrt(size(b, 2)) * delta;
        cvx_end
end
end




