function [epind, epindlong] = find_vicinity(sdl, curryloc, currytri, method, param)

currytri(:, sum(currytri > size(curryloc, 2)) ~= 0) = [];
[~, sdind]  = find_nvoxel(sdl, curryloc); nnd  = length(sdind);
epind  = cell(nnd, 1); % entire patch indices
epa  = zeros(1, nnd);


for i_nnd = 1 : nnd
    
    %=====================================================================%
    switch method
        
        %=================================================================%
        case 'neighbourhood'
            epind{i_nnd} = sdind(i_nnd);
            for i_prm = 1 : param(i_nnd)
                tri = currytri(:, find_triind(currytri,epind{i_nnd},1:3));
                epind{i_nnd}  = unique(tri(:).');
            end
            id3  = find_triind(currytri, epind{i_nnd}, 3);
            epa(i_nnd) = find_area(curryloc, currytri(:, id3));
            
        %=================================================================%
        case 'area'
            epind{i_nnd} = sdind(i_nnd);
            oldind  = [];
            while( epa(i_nnd) < param(i_nnd) )
                tri = currytri(:, find_triind(currytri,epind{i_nnd},1:2));
                for i_tri = 1 : size(tri, 2)
                    tempind = unique([epind{i_nnd} tri(:, i_tri).']);
                    if(isequal(tempind, epind{i_nnd})), continue; else
                        epind{i_nnd}  = tempind; 
                        newind  = find_triind(currytri, epind{i_nnd}, 3);
                        newtri  = currytri(:, setdiff(newind, oldind));
                        oldind  = newind;
                        epa(i_nnd) = epa(i_nnd)+find_area(curryloc,newtri);
                        if( epa(i_nnd) >= param(i_nnd) ), break; end
                    end
                end
            end

        %=================================================================%
        otherwise
            display('================== Weird Request ==================')
        
    end
    
    %=====================================================================%
    epind(i_nnd)  = revise_sdp(epind(i_nnd), sdind(i_nnd));
    
end

epindlong = ppatches(epind, 'row');

end
