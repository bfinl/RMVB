function [KLD] = anlys_KLD(currycdr, curryloc, source, Phi)

bglevel  = 0.01; % background level
[P, T]   = findpeaks(std(Phi, 1, 1)); 
T        = T(P >= max(P)/sqrt(10));

%%
for t = 1 : length(T)
    cdr  = norms(squeeze(currycdr(4, :, :)), 2, 2).';
%     cdr  = std(squeeze(currycdr(4, :, :)), 1, 2).';
%     cdr  = abs(squeeze(currycdr(4, :, T(t))));

    if(sum(cdr ~= 0))
        cdr(cdr == 0)  = bglevel * min(cdr(cdr ~= 0)) / sum(cdr == 0);

        %%
        ref  = zeros(size(cdr));
        for i = 1 : length(source.epl)
            [~, ind]  = find_nvoxel(source.epl{i}, curryloc); 
            ref(ind)  = norms(source.sdm(i, :), 2, 2);
%             ref(ind)  = std(source.sdm(i, :), 1, 2);
%             ref(ind)  = abs(source.sdm(i, T(t)));
        end
        ref(ref == 0)  = bglevel * min(ref(ref ~= 0)) / sum(ref == 0);

        %%
        KLD(t)  = (KLDiv(ref, cdr)) / 2;
        % KLD  = (KLDiv(ref, cdr) + KLDiv(cdr, ref)) / 2;

    else 
        KLD(t)  = inf;
    end
end

KLD  = mean(KLD);

end


function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
% downloaded from Mathworks File Exchange


if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end

end
