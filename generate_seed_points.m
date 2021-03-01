function [SDI, TAG] = generate_seed_points(Phi, LFD, ULFD, TRI, OND, NPL, sdp_type, lfd_type, varargin) 
% SDL : seed location, SDI : seed index
% NHD : neighbourhood, NPL : number of parcels
if(strcmp(lfd_type, 'rot')), d = 3; elseif(strcmp(lfd_type, 'fxd')), d = 1; end

%%
switch sdp_type
    case 'GROVA'
        Phi  = Phi ./ repmat(norms(Phi), [size(Phi, 1), 1]);
        LFD  = LFD ./ repmat(norms(LFD), [size(LFD, 1), 1]);
        [U, ~, ~]  = svds(Phi, rank(Phi));
        MSP  = abs(LFD.'*U).^2; 
        MSP  = reshape(MSP, d, []); 
        MSP  = sum(MSP, 1);
        MSP  = reshape(MSP, size(LFD, 2)/d, []); 
        [m, n]  = size(MSP);
        [MSPMAX, TAG]  = max(MSP, [], 2);
        MSP(MSP < repmat(MSPMAX, [1, n])) = -1;

    case 'MUSIC'
        [U, ~, ~]  = svds(Phi, rank(Phi));  
        MSP  = zeros(size(LFD, 2)/3, 1);
        count  = 0;
        for i = 1 : d : size(LFD, 2)
            count  = count + 1;
            MSP(count)  = max(eig((ULFD(:, i:i+d-1).') * U * (U.') * ULFD(:, i:i+d-1)));
        end
        [m, n]  = size(MSP);
        TAG  = ones(m, 1);
        
    case 'OFFSP'
        [U, ~, ~]  = svds(Phi, rank(Phi));               
        MSP  = zeros(size(LFD, 2)/3, size(U, 2));
        for i_U = 1 : size(U, 2)
            count  = 0;
            for i = 1 : d : size(LFD, 2)
                count  = count + 1;
                MSP(count, i_U)  = max(eig((ULFD(:, i:i+d-1).') * U(:, i_U) * (U(:, i_U).') * ULFD(:, i:i+d-1)));
            end
        end
        [m, n]  = size(MSP);
        TAG  = ones(m, 1);
        
    case 'EXTNL'
        MSP  = varargin{:};
        [m, n]  = size(MSP);
        TAG  = ones(m,  1);
        
end

%%
SDI  = [];
nclmn  = 1;
while(length(SDI)<NPL && ~isequal(MSP, -1*ones([m, n])))
    [~, idx] = max(MSP(:, nclmn));
    SDI = [SDI, idx];

    for i_nhd = 1 : OND
        tri  = TRI(:, find_triind(TRI, idx, 1:3));
        idx  = unique(tri(:));
    end
    MSP(ismember(1:m, idx), nclmn) = -1;

    if(isequal(MSP(:, nclmn), -1*ones([m, 1])))
        nclmn  = nclmn + 1;
    end
end

%%
if(length(SDI) ~= NPL)
    fprintf('Warning: %f out of %f parcels were reached!\n', length(SDI), NPL)
end

end