function ind  =  extend_ind(ind, varargin)

if(isempty(varargin) || strcmp(varargin{:}, 'rot'))
    offset  = repmat([0; 1; 2], [1, length(ind)]);
    ind  = 3 * ind -2;
    ind  = repmat(ind, [3, 1]) + offset;
    ind  = ind(:);
end

end