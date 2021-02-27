function [EPI, PMP, OND] = generate_parcellation(SDI, TAG, TRI, M)
% SDI: seed index
% TAG: TAG for each dipole
% TRI: tri matrix
%   N: number of valid dipoles
% EOI: entire patch indices
% PMP: parcellation map
% OND: oredr of neighbourhood


%%
N = max(unique(TRI(:)));
EPI = cell(length(SDI), 1);
for i = 1 : length(SDI), EPI{i} = SDI(i); end
DAS = ones(1, N); % dipole availability status
for i_epi = 1 : length(EPI), DAS(EPI{i_epi}) = -1; end


%%
count =0;
OND = zeros(length(EPI), 1);
while(~isequal(DAS, -1*ones(1, N)))
    for i_epi = 1 : length(EPI)
        ind  = find_triind(TRI, EPI{i_epi}, 1:3);
        tri  = TRI(:, ind); tri  = tri(:); 
%         tri  = tri(TAG(tri) == TAG(SDI{i_epi}));
        tri  = tri(DAS(tri) ~= -1);

        if(~isempty(tri)), OND(i_epi)=OND(i_epi)+1; end
        EPI{i_epi}  = unique([EPI{i_epi}, tri.']);
        DAS(EPI{i_epi})  = -1;
        TRI(:, ind) = [];
    end
    count = count+1
end
EPI = revise_sdp(EPI, SDI);


%% removes the off-inner skull indices
for i_epi = 1 : length(EPI)
    EPI{i_epi}(EPI{i_epi} > M)  = [];
end

    
%%
PMP  = zeros(M, 1);
CLR  = randperm(length(EPI), length(EPI)); 
for i = 1 : length(EPI), PMP(EPI{i}) = CLR(i); end


%%
epi = ppatches(EPI, 'row');
len = 0; for i = 1 : length(EPI), len  = len + length(EPI{i}); end
if(abs(len-M) || abs(length(epi)-M)), disp('Warning: sth went wrong!'), end
fprintf('Average region growing iterations: %f\n', mean(OND));


end