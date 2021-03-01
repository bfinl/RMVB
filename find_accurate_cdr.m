function [mnt]  = find_accurate_cdr(sdl, epl, Phi_bl, Phi_az, curryloc, currylfd, lfd_type, cdr_type, dsply_opt, frq, dli, varargin)


switch cdr_type
    
    case 'std' % standard
        currycdr  = varargin{:};
        [~, ind]  = find_nvoxel(sdl, curryloc);
        mnt(1:length(ind), :) = squeeze(currycdr(4, ind, :));
    
    case 'ext' % extent analysis
        currycdr  = varargin{:};
        mnt  = zeros(length(epl), size(currycdr, 3));
        for i = 1 : length(epl)
            [~, ind]  = find_nvoxel(epl{i}, curryloc);
            tmprycdr  = revise_cdr_ori(currycdr(:, ind, :), 'max');
            mnt(i, :) = squeeze(mean(tmprycdr(4, :, :), 2));
        end

    otherwise % solves the problem again using some algorithm
        [~, ind]  = find_nvoxel(sdl, curryloc);
        currylfd  = currylfd(:, extend_ind(ind, lfd_type)); 
        currycdr  = solve_inverse_problem(Phi_bl, Phi_az, currylfd, lfd_type, cdr_type);    
        tmprycdr  = revise_cdr_ori(currycdr, 'max');
        mnt(1:length(ind), :) = squeeze(tmprycdr(4, :, :));

end

      

%%

if( strcmp(dsply_opt, 'Display') )
    nnd  = size(sdl, 2);
    nts  = size(Phi_az, 2);
    for i = 1 : nnd
        lgnd{i} = ['Source ' num2str(dli(i))]; 
    end
    tp = 0:1/frq:(nts-1)/frq; tp = repmat(tp.', 1, nnd);
    h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
    plot(tp, mnt.'), legend(lgnd), xlabel('Time (ms)'), ylabel('Moment (\muAmm)')
end


end






