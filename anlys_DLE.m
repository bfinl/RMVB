function [est_sdl, est_epl, DLE, DLI] = anlys_DLE(curryloc, currycdr, currytri, source, th, th_olp, nnd_type, loc_type)

r  = 12.5;
est_sdl  = [];
est_epl  = {};
DLE  = [];
DLI  = [];

%% Number of nodes estimation

switch nnd_type
   
    case 'Predefined'
        
        est_nnd  = source.nnd;
        
        
    case 'PCA'
        
        % ...
        
end


%% Main (localization)

switch loc_type
    
    case 'MAX-FOC'
        amp  = abs(squeeze(currycdr(4, :, :))).^1;
        amp  = max(amp, [], 2);
        th_ind  = find(amp > 0 * max(amp));
        [est_val, ~, ind]  = find_local_max(amp(th_ind), curryloc(:, th_ind), r);
        est_ind  = th_ind(ind);
        [est_sdl, est_epl] = find_entities(est_val, est_ind, est_nnd, currycdr, curryloc, currytri, th_olp, 'FDA');
        tru_sdl  = source.sdl;
     
    case 'VAR-FOC'
        amp  = squeeze(currycdr(4, :, :));
%         amp  = revise_cdr_dir(currycdr);
%         amp  = abs(amp(4, :, :));
        amp  = var(amp, 1, 2);
        th_ind  = find(amp > 0 * max(amp));
        [est_val, ~, ind]  = find_local_max(amp(th_ind), curryloc(:, th_ind), r);
        est_ind  = th_ind(ind);
        [est_sdl, est_epl] = find_entities(est_val, est_ind, est_nnd, currycdr, curryloc, currytri, th_olp, 'FDA');
        tru_sdl  = source.sdl;
        
    case 'PCA-FOC'
        est_ind  = [];
        est_val  = [];
        J  = squeeze(currycdr(4, :, :));
        [U, ~, ~]  = svds(J, min(est_nnd, rank(J)));
        for i = 1 : size(U, 2)
            amp = abs(U(:, i));
            th_ind  = find(amp > th * max(amp));
            [val, ~, ind]  = find_local_max(amp(th_ind), curryloc(:, th_ind), r);
            est_ind  = [est_ind; th_ind(ind)]; 
            est_val  = [est_val; val];
        end
        [~, ~, ind]  = find_local_max(est_val, curryloc(:, est_ind), r);
        est_val  = est_val(ind);
        est_ind  = est_ind(ind);
        [est_sdl, est_epl] = find_entities(est_val, est_ind, est_nnd, currycdr, curryloc, currytri, th_olp, 'FDA');
        tru_sdl  = source.sdl;
 
    case 'PCA-EXT'
        est_ind  = [];
        est_val  = [];
        J  = squeeze(currycdr(4, :, :));
        [U, ~, ~]  = svds(J, min(est_nnd, rank(J)));
        for i = 1 : size(U, 2)
            amp = abs(U(:, i));
            th_ind  = find(amp > th * max(amp));
            [val, ~, ind]  = find_local_max(amp(th_ind), curryloc(:, th_ind), r);
            est_ind  = [est_ind; th_ind(ind)]; 
            est_val  = [est_val; val];
        end
        [~, ~, ind]  = find_local_max(est_val, curryloc(:, est_ind), r);
        est_val  = est_val(ind);
        est_ind  = est_ind(ind);
        [est_sdl, est_epl, ~, est_epi] = find_entities(est_val, est_ind, est_nnd, currycdr, curryloc, currytri, th_olp, 'EDA');
%         est_sdl  = find_center(est_epl, est_epi, currycdr);
        est_sdl  = find_center(est_epl);
        tru_sdl  = find_center(source.epl);

end



%% dipole localization error

est_nnd   = size(est_sdl, 2);
dist_mat  = zeros(est_nnd, source.nnd);        
for i = 1 : est_nnd
    for j = 1 : source.nnd
        dist_mat(i, j) = norm(est_sdl(:, i) - tru_sdl(:, j));
    end 
end
[rule, cost]  = assignmentoptimal(dist_mat);
% ave_err  = cost / min(est_nnd, source.nnd);

rule(rule == 0) = [];
[DLI, ind] = sort(rule, 'ascend');
est_sdl  = est_sdl(:, ind);
est_epl  = est_epl(1, ind);


for i = 1 : length(ind)
    DLE(i)  = norm(est_sdl(:, i) - tru_sdl(:, DLI(i)));
end


end


function [est_sdl, est_epl, est_sdi, est_epi] = find_entities(est_val, est_ind, est_nnd, currycdr, curryloc, currytri, th_olp, type)

[~, IDX]  = sort(est_val, 'descend');
if( length(IDX) < est_nnd )
    fprintf('WARNING: %d node(s) out of %d recovered! \n', length(IDX), est_nnd)
end
est_nnd  = min(length(IDX), est_nnd);
est_sdi  = est_ind(IDX(1:est_nnd));
est_sdl  = curryloc(:, est_sdi);
est_epi  = cell(1, length(est_sdi));
est_epl  = cell(1, length(est_sdi));
for i = 1 : length(est_sdi)
    est_epi{1, i}  = est_sdi(i);
    est_epl{1, i}  = est_sdl(:, i);
end
        
switch type
    
    case 'FDA' % focal dipole analysis
        % no need to do anything
        
    case 'EDA' % extent dipole analysis
        cdr  = norms(squeeze(currycdr(4, :, :)), 2, 2);
        cdr  = [cdr; zeros(max(currytri(:))-size(curryloc, 2), 1)];
        das  = ones(1, length(cdr)); % dipole availability status
        das(cdr < th_olp*max(cdr))  = 0; 
        
        % main
        BLEAN  = 1;
        while(BLEAN)
            BLEAN  = 0;
            for ii = 1 : length(est_sdi)  
                tri  = currytri(:, find_triind(currytri, est_epi{ii}, 1:2));
                ind  = unique(tri(:)).';
                ind  = setdiff(ind, est_epi{ii});
                ind  = ind(das(ind) ~= 0);
                est_epi{ii}  = [est_epi{ii} ind];
                BLEAN  = BLEAN | ~isempty(ind);
                das(est_epi{ii}) = 0;
            end
        end
        est_epi = revise_sdp(est_epi, est_sdi);
        
        % removes the off-inner skull indices
        for i_epi = 1 : length(est_epi)
            est_epi{i_epi}(est_epi{i_epi} > size(curryloc, 2))  = [];
            est_epl{1, i_epi}  = curryloc(:, est_epi{1, i_epi});
        end
        
        
end

end


function ctr_loc = find_center(epl, varargin)

if(isempty(varargin))
    
    ctr_loc  = zeros(3, length(epl));
    for i = 1 : length(epl)
        ctr_loc(:, i)  = find_nvoxel(mean(epl{i}, 2), epl{i});
    end
    
else
    
    epi  = varargin{1};
    cdr  = varargin{2};
    www  = norms(squeeze(cdr(4, :, :)), 2, 2);
    
    ctr_loc  = zeros(3, length(epl));
    for i = 1 : length(epl)
        ctr_out = epl{i}*www(epi{i})/sum(www(epi{i}));
        ctr_loc(:, i)  = find_nvoxel(ctr_out, epl{i});
    end
    
end

end
