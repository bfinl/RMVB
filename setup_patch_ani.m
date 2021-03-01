function [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch_ani(curryloc, curryori, currytri, curryani, nnd, sdm)

epl  = cell(1, nnd);
epo  = cell(1, nnd);
epm  = cell(1, nnd);
epa  = zeros(1, nnd);
sdl  = zeros(3, nnd);
sdo  = zeros(3, nnd);
sdi  = randperm(length(curryani), nnd); 

for i_nnd = 1 : nnd
    
    sdl(:, i_nnd)  = curryloc(:, curryani{sdi(i_nnd)}(1));
    sdo(:, i_nnd)  = curryori(:, curryani{sdi(i_nnd)}(1));
    epl{1, i_nnd}  = curryloc(:, curryani{sdi(i_nnd)});
    epo{1, i_nnd}  = curryori(:, curryani{sdi(i_nnd)});
    
    epm{1, i_nnd}  = repmat(sdm(i_nnd, :), [length(curryani{sdi(i_nnd)}), 1]); 
    ind  = find_triind(currytri, curryani{sdi(i_nnd)}, 3);
    epa(1, i_nnd) = find_area(curryloc, currytri(:, ind));

end

epc  = [];
% tmp  = ppatches(epm, 'col'); tim  = size(tmp, 2);
% tmp  = tmp - repmat(mean(tmp, 2), [1, tim]);
% epc  = (tmp * tmp.') / tim; % covariance

end