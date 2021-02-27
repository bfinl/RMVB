function [LoMaVal, LoMaLoc, LoMaInd] = find_local_max(J, curryloc, r)

LoMaInd = [];
m = length(J);

continue_flag  = ones(1, m);
for i = 1 : m
    if( continue_flag(i) )
        dist = norms(curryloc - repmat(curryloc(:, i), [1, m]));
        adjVec = (dist <= r);
        blean = nnz((J(i) - J(adjVec)) < 0);
        if( ~blean )
            LoMaInd = [LoMaInd i];
            continue_flag(adjVec)  = 0;
        end
    end
end

LoMaVal = J(LoMaInd);
LoMaLoc = curryloc(:, LoMaInd);

end


% for i = 1 : m
%     dist = norms(curryloc - repmat(curryloc(:, i), [1, m]));
%     adjVec = (dist <= r) & (dist > 0);
%     blean = nnz((J(i) - J(adjVec)) <= 0);
%     blean_curry  = ismember(curryloc(:, i).', curryloc(:, LoMaInd).', 'rows'); 
%     if( ~blean && ~blean_curry )
%         LoMaInd = [LoMaInd i];
%     end
% end


% continue_flag  = ones(1, m);
% [~, ind_srt]  = sort(J, 'descend');
% for i = 1 : m
%     if( continue_flag(i) )
%         dist = norms(curryloc - repmat(curryloc(:, ind_srt(i)), [1, m]));
%         adjVec = (dist <= r);
%         blean = nnz((J(ind_srt(i)) - J(adjVec)) < 0);
%         if( ~blean )
%             LoMaInd = [LoMaInd ind_srt(i)];
%             continue_flag((dist(ind_srt) <= r))  = 0;
%         end
%     end
% end