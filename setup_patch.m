function [sdl, sdo, epl, epo, epm, epa, epc] = setup_patch(curryloc, curryori, currytri, nnd, sdm, method, param, varargin)

epl  = cell(1, nnd);
epo  = cell(1, nnd);
epm  = cell(1, nnd);
epa  = zeros(1, nnd);
sdl  = zeros(3, nnd);
sdo  = zeros(3, nnd);
epind  = cell(nnd, 1); % entire patch indices
sdind  = zeros(nnd, 1); % seed indices
default_area_range  = pi * (10:30) .^ 2; 
currytri(:, sum(currytri > size(curryloc, 2)) ~= 0) = [];
if(isempty(varargin)), refst  = 1:size(curryloc, 2); else refst = varargin{:}; end 
if(strcmp(param, 'random')), param  = default_area_range(randperm(length(default_area_range), nnd)); end

for i_nnd = 1 : nnd
    
    %=====================================================================%
    switch method
        
        %=================================================================%
        case 'neighbourhood'
            sdind(i_nnd) = refst(randi(length(refst))); 
            epind{i_nnd} = sdind(i_nnd);
            for i_prm = 1 : param(i_nnd)
                tri = currytri(:, find_triind(currytri,epind{i_nnd},1:3));
                epind{i_nnd}  = unique(tri(:));
            end
            id3  = find_triind(currytri, epind{i_nnd}, 3);
            epa(i_nnd) = find_area(curryloc, currytri(:, id3));
            
        %=================================================================%
        case 'area'
            sdind(i_nnd) = refst(randi(length(refst))); 
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
        case 'area_slow'
            sdind(i_nnd) = refst(randi(length(refst))); 
            epind{i_nnd} = sdind(i_nnd);           
            oldind  = [];
            while( epa(i_nnd) < param(i_nnd) )
                tri = currytri(:, find_triind(currytri,epind{i_nnd},1:2));
                for i_tri = 1 : size(tri, 2)
                    tempind = unique([epind{i_nnd} tri(:, i_tri).']);
                    if(isequal(tempind, epind{i_nnd})), continue; else
                        epind{i_nnd}  = tempind; 
                        ind  = find_triind(currytri, epind{i_nnd}, 3);
                        tri  = currytri(:, ind);
                        epa(i_nnd) = find_area(curryloc, tri);
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
    refst  = setdiff(refst, epind{i_nnd});
    sdl(:, i_nnd)  = curryloc(:, sdind(i_nnd));
    sdo(:, i_nnd)  = curryori(:, sdind(i_nnd)); 
    epl{1, i_nnd}  = curryloc(:, epind{i_nnd});
    epo{1, i_nnd}  = curryori(:, epind{i_nnd}); 
    epm{1, i_nnd}  = repmat(sdm(i_nnd, :), [length(epind{i_nnd}), 1]); 
    
end

epc  = find_cov(ppatches(epm, 'col'), 'SC');

end
