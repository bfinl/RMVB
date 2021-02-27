function [C, lglik] = find_cov(X, type, varargin)

if(isempty(varargin))
    remove_mean_opt = 'remove mean';
else 
    remove_mean_opt = varargin{:};
end

T  = size(X, 2);
if( strcmp(remove_mean_opt, 'remove mean') )
    X  = X - repmat(mean(X, 2), [1, T]);
end
S  = (X * X.') / T;

switch type
    
    %*********************************************************************%
    case 'EM' % Empirical
        C  = S;

    %*********************************************************************%
    case 'DG' % Diagonal 
        C  = diag(diag(S));
        
    %*********************************************************************%
    case 'SC' % Shrunk Covariance
        kfold  = 3;
        alpha  = logspace(log10(1e-2), log10(1e0), 30);
        cvprtn = cvpartition(T, 'KFold', kfold);
        lglik  = zeros(length(alpha), kfold); 
        for i = 1 : length(alpha)
            for j = 1 : kfold
                cvind  = training(cvprtn, j);
                i_trn  = cvind == 1;
                i_vld  = cvind == 0;
                S_trn  = (X(:, i_trn) * X(:, i_trn).') / length(i_trn);
                M_trn  = mean(diag(S_trn));
                C_trn  = (1-alpha(i))*S_trn + alpha(i)*M_trn*eye(size(X, 1));
                S_vld  = (X(:, i_vld) * X(:, i_vld).') / length(i_vld);
                lglik(i, j)  = - trace(S_vld * (C_trn^-1)) - sum(log(eig(C_trn))); 
            end
        end
        [~, ind]  = max(mean(lglik, 2));
        alpha_cv  = alpha(ind);
        M  = mean(diag(S));
        C  = (1-alpha_cv)*S + alpha_cv*M*eye(size(X, 1));

    %*********************************************************************%
    case 'FA' % Factor Analysis        
        kfold  = 3;
        erank  = 1:5:128;  
        cvprtn = cvpartition(T, 'KFold', kfold);
        lglik  = zeros(length(erank), kfold); 
        for i = 1 : length(erank)
            for j = 1 : kfold
                cvind  = training(cvprtn, j);
                i_trn  = cvind == 1;
                i_vld  = cvind == 0;
                [F, psi] = FA_abbas(X(:, i_trn), erank(i), 100, 1e-3);
                C_trn  = F*(F.') + diag(psi);
                S_vld  = (X(:, i_vld) * X(:, i_vld).') / length(i_vld);
                lglik(i, j)  = - trace(S_vld * (C_trn^-1)) - sum(log(eig(C_trn))); 
            end
        end
        [~, ind]  = max(mean(lglik, 2));
        rank_cv   = erank(ind);
        [F, psi]  = FA_abbas(X, rank_cv, 100, 1e-3);
        C  = F*(F.') + diag(psi);       
        
    %*********************************************************************%
    case 'CD' % Combined
        lglik  = zeros(1, 4);
        C  = zeros(size(X, 1), size(X, 1), 4); 
        [C(:, :, 1), lglik(1)] = find_cov(X, 'EM');
        [C(:, :, 2), lglik(2)] = find_cov(X, 'DG');
        [C(:, :, 3), lglik(3)] = find_cov(X, 'SC');
        [C(:, :, 4), lglik(4)] = find_cov(X, 'FA');
        [~, ind]  = max(lglik);
        C  = C(:, :, ind);

end

lglik  = - trace(S*(C^-1)) - sum(log(eig(C))); 

end