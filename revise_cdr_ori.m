function [J, ave_ori]  = revise_cdr_ori(currycdr, revision_type, varargin)

switch revision_type
    
    case 'ext'
        ave_ori  = varargin{:};
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* sign(dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
    case 'max'        
        [~, ind]  = max(squeeze(abs(currycdr(4, : , :))), [], 2);
        for i = 1 : size(currycdr, 2)
            ave_ori(:, i)  = currycdr(1:3, i, ind(i));
        end
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* sign(dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
    case 'svd-v1'
        ave_ori  = zeros(3, size(currycdr, 2));
        for i = 1 : size(currycdr, 2)
            [U, S, ~]  = svds(squeeze(currycdr(1:3, i , :)), 3);
            ave_ori(1:3, i)  = U(:, 1);
        end
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* sign(dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
    case 'svd-v2'
        ave_ori  = zeros(3, size(currycdr, 2));
        for i = 1 : size(currycdr, 2)
            [U, S, ~]  = svds(squeeze(currycdr(1:3, i , :)), 3);
            ave_ori(1:3, i)  = U(:, 1);
        end
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* (dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
    case 'mean'
        pow  = varargin{:};
        ave_ori  = mean(repmat(currycdr(4, : , :).^pow, [3, 1, 1]) .* currycdr(1:3, : , :), 3);
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* sign(dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
    case 'median'
        ave_ori  = median(currycdr(1:3, : , :), 3);
        ave_ori  = ave_ori ./ repmat(norms(ave_ori), [3, 1]);
        J(1:3, :, :)  = repmat(ave_ori, [1, 1, size(currycdr, 3)]);
        J(4, :, :)  = currycdr(4, : , :) .* sign(dot(J(1:3, :, :),  currycdr(1:3, :, :), 1));
        J(isnan(J)) = 0;
        
end

end