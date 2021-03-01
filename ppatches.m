function patchedpatches = ppatches(patches, type)
% this function patches the patches :)

patchedpatches  = [];
switch type
    case 'row'
        for i = 1 : length(patches)
            patchedpatches  = [patchedpatches patches{i}];
        end
        
    case 'col'
        for i = 1 : length(patches)
            patchedpatches  = [patchedpatches; patches{i}];
        end
end

end