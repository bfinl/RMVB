function K  = revise_lfd_dir(currylfd, varargin)

if(isempty(varargin))
    lfd_typ  = 'rot';
else
    lfd_typ  = varargin{:};
end

switch lfd_typ
    
    case 'rot'
        K  = repmat(diag(currylfd(4:6, 1:3)).', size(currylfd, 1)-6, size(currylfd, 2)/3) .* currylfd(7:end, :); 
        
    case 'fxd'
        K  = currylfd(7:end, :);
end
        

end