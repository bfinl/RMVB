function boundary_edges = find_boundary_edges(rsctntri)
edge_map = [];
for i = 1 : size(rsctntri, 2)
    for j = 1 : 3
        edge_map  = [edge_map [rsctntri(j, i); rsctntri(mod(j, 3)+1, i)]];
    end
end

rep_num  = zeros(1, size(edge_map, 2));
for i = 1 : size(edge_map, 2)
   rep_num(i)  = nnz(ismember(edge_map.', edge_map(:, i).', 'rows') | ...
                     ismember(edge_map.', edge_map(2:-1:1, i).', 'rows')); 
end

boundary_edges  = edge_map(:, rep_num == 1); % boundaries are the edges belonging to only one facet 



end