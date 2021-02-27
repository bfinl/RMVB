function [bin_ctr, bin_val, bin_err] = my_binstats(dis_rec, err_rec, nbins, type)

mn  = min(dis_rec);
mx  = max(dis_rec);
bin_edg  = linspace(mn, mx, nbins+1);
bin_ctr  = zeros(nbins, 1);
bin_val  = zeros(nbins, 1);
bin_err  = zeros(nbins, 2);

switch type
        
    case 'mean'
        for i = 1 : nbins
            ind  = find(bin_edg(i)<=dis_rec & dis_rec<=bin_edg(i+1));
            bin_ctr(i)  = (bin_edg(i) + bin_edg(i+1))/2;
            bin_val(i)  = mean(err_rec(ind));
            bin_err(i, 1)  = 0;
            bin_err(i, 2)  = std(err_rec(ind), 1);
        end
    
    case 'median'
        for i = 1 : nbins
            ind  = find(bin_edg(i)<=dis_rec & dis_rec<=bin_edg(i+1));
            bin_ctr(i)  = (bin_edg(i) + bin_edg(i+1))/2;
            dataquntl   = quantile(err_rec(ind), [0.25 0.5 0.75]);
            bin_val(i)  = dataquntl(2);
            bin_err(i, 1)  = dataquntl(2) - dataquntl(1);
            bin_err(i, 2)  = dataquntl(3) - dataquntl(2);
        end
        
end

end