function triind = find_triind(tri, locind, m, varargin)

typ  = 'regular';
ind  = ismember(tri(:), locind);
ind  = reshape(ind, size(tri, 1), []);
if(~isempty(varargin)), typ = varargin{:}; end;

switch typ
    case 'regular'
        triind  = zeros(1, size(tri, 2));
        for i = length(m) : -1 : 1
            triind  = triind | (sum(ind) == m(i));
        end
        triind  = find(triind == 1);

    case 'priority'
        triind  = [];
        for i = length(m) : -1 : 1
            triind  = [triind find(sum(ind) == m(i))];
        end
end

end