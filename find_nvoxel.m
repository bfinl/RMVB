function [locVal, locInd] = find_nvoxel(loc_mat, curryloc)

locVal  = [];
locInd  = [];

for i = 1 : size(loc_mat, 2)
    
    loc            =  loc_mat(:, i); 
    loc            =  repmat(loc, 1, size(curryloc, 2));
    [~, idx]       =  min(norms(curryloc-loc));
    locInd         =  [locInd idx];
    locVal         =  [locVal curryloc(:, idx)];

end


end