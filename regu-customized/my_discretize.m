function [bin_ctr, bin_ave, bin_err] = my_discretize(dis_rec, err_rec, fitted_data, nbins, type)

mn  = min(dis_rec);
mx  = max(dis_rec);
bin_edg  = linspace(mn, mx, nbins+1);
bin_ctr  = zeros(nbins, 1);
bin_ave  = zeros(nbins, 1);
bin_err  = zeros(nbins, 2);

switch type
    
    case 'dscrtz-v1'
        for i = 1 : nbins
            ind  = find(bin_edg(i)<=dis_rec & dis_rec<=bin_edg(i+1));
            err_fit  = err_rec(ind).' - fitted_data(dis_rec(ind));
            err_mns  = err_fit(err_fit <= 0);
            err_pls  = err_fit(err_fit >= 0);
            bin_ctr(i)  = (bin_edg(i) + bin_edg(i+1))/2;
            bin_ave(i)  = fitted_data(bin_ctr(i));
            bin_err(i, 1)  = max(0, bin_ctr(i)-norm(err_mns) / sqrt(length(err_mns)));
            bin_err(i, 2)  = norm(err_pls) / sqrt(length(err_pls));
        end
        
    case 'dscrtz-v2'
        for i = 1 : nbins
            ind  = find(bin_edg(i)<=dis_rec & dis_rec<=bin_edg(i+1));
            bin_ctr(i)  = (bin_edg(i) + bin_edg(i+1))/2;
            bin_ave(i)  = mean(err_rec(ind));
            bin_err(i, 1)  = std(err_rec(ind), 1);
            bin_err(i, 2)  = bin_err(i, 1);
        end
        
end

end