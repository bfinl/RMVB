function orthglfd = find_orthg_lfd(currylfd, varargin)

if(isempty(varargin))
    lfd_typ  = 'rot';
else
    lfd_typ  = varargin{:};
end

switch lfd_typ
    
    case 'rot'
        for i = 1 : 3 : size(currylfd, 2)
            [U_k, ~, ~]  = svds(currylfd(:, i:i+2), 3);
            orthglfd(:, i:i+2)  = U_k;
        end
        
    case 'fxd'
        orthglfd  = currylfd ./ repmat(norms(currylfd), [size(currylfd, 1), 1]);
        
end

end