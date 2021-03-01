function area = find_area(curryloc, tri)

area  = 0;

for i = 1 : size(tri, 2)
    acor  = curryloc(:, tri(1, i));
    bcor  = curryloc(:, tri(2, i));
    ccor  = curryloc(:, tri(3, i));
    
    asid  = norms(acor-bcor);
    bsid  = norms(bcor-ccor);
    csid  = norms(ccor-acor);
    
    hper  = (asid + bsid + csid) / 2;
    area  = area + sqrt(hper * (hper-asid) * (hper-bsid) * (hper-csid));
end

end



% function area = find_area(curryloc, tri)
% 
% aloc  = curryloc(:, tri(1, :));
% bloc  = curryloc(:, tri(2, :));
% cloc  = curryloc(:, tri(3, :));
% area  = zeros(1, size(tri, 2));
% 
% parfor i = 1 : size(tri, 2)
% 
%     acor  = aloc(:, i);
%     bcor  = bloc(:, i);
%     ccor  = cloc(:, i);
%     
%     asid  = norms(acor-bcor);
%     bsid  = norms(bcor-ccor);
%     csid  = norms(ccor-acor);
%     
%     hper  = (asid + bsid + csid) / 2;
%     area(i)  = sqrt(hper * (hper-asid) * (hper-bsid) * (hper-csid));
% end
% area  = sum(area);
% 
% end