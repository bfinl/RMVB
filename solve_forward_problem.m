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

function Phi = solve_forward_problem(currylfd, curryloc, source, lfd_type )
% This function solves the forward problem to generate electrod recordings 
% corresponding to a given source activity in the brain. 
% 
% Parmeters:
%     currylfd: leadfield matrix
%     curryloc: 3d locations of the brain vertices
%     source: underlying source activity
%     lfd_type: either 'rot' or 'fxd' for a rotational or fixed-orientation
%     model
% 
% Returns:
%     Phi: The generated electrode recordings
%         
% Written by Seyed Amir Hossein Hosseini
% 09/01/2017


switch lfd_type 
    
    case 'rot'
        epm = [];
        epl = ppatches(source.epl, 'row');
        epo = ppatches(source.epo, 'row');
        for i = 1 : length(source.epl)
            epm  = [epm; repmat(source.sdm(i, :), [size(source.epl{i}, 2), 1])];
        end
        [~, ind]  = find_nvoxel(epl, curryloc);
        ind = extend_ind(ind);

        Phi = [];
        for i = 1 : source.nts
            J_i  = (epo .* repmat(epm(:, i).', [3, 1])); J_i  = J_i(:);
            Phi  = [Phi currylfd(:, ind)*J_i];
        end

        
    case 'fxd'
        
        epm = [];
        epl = ppatches(source.epl, 'row');
        for i = 1 : length(source.epl)
            epm  = [epm; repmat(source.sdm(i, :), [size(source.epl{i}, 2), 1])];
        end
        [~, ind]  = find_nvoxel(epl, curryloc);
        Phi = currylfd(:, ind) * epm;        
        
        
    otherwise
        display('*************** Weird Request ***************')
        

end