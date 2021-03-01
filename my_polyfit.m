function [fitted_data, ind_mdl] = my_polyfit(xdata, ydata, nmodels, kfold, varargin)

if(isempty(varargin))
    
    cvprtn  = cvpartition(length(xdata), 'KFold', kfold);
    mdlflderr  = zeros(nmodels, kfold); 
    for i_mdl = 1 : nmodels
        for i_fld = 1 : kfold
            cvind  = training(cvprtn, i_fld);
            ind_trn  = cvind == 1;
            ind_vld  = cvind == 0;
            LB  = -inf * ones(i_mdl+1, 1); LB(end) = 0;
            UB  = +inf * ones(i_mdl+1, 1); UB(end) = 0;
            fitted_data  = fit(xdata(ind_trn), ydata(ind_trn), ['Poly' num2str(i_mdl)], 'Lower', LB, 'Upper', UB);
            mdlflderr(i_mdl, i_fld)  = mean((ydata(ind_vld) - fitted_data(xdata(ind_vld))).^2);
        end
    end
    [~, ind_mdl]  = min(mean(mdlflderr, 2));
    
else
    
    ind_mdl  = varargin{:};
    
end

LB  = -inf * ones(ind_mdl+1, 1); LB(end) = 0;
UB  =  inf * ones(ind_mdl+1, 1); UB(end) = 0;
fitted_data  = fit(xdata, ydata, ['Poly' num2str(ind_mdl)], 'Lower', LB, 'Upper', UB, 'Robust','on');


end