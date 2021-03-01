function [gamma2, rel_err] = anlys_CNN(ts, pMin, pMax, lFreq, hFreq, Fs, shuffTimes, sigLevel, nlz_opt, arfit_opt, dsply_opt, varargin)


nnd  = size(ts, 1);
nts  = size(ts, 2);



%%
if (strcmp(nlz_opt, 'Normalize')) 
    ts  = (ts - repmat(mean(ts, 2), [1, nts])) ./ repmat(std(ts, 1, 2), [1, nts]);
end
ts = ts.';



%%
[~, ~, ~, sbc, fpe, ~] =    arfit(ts, pMin, pMax, arfit_opt);         
[~, ind_sbc]           =    min(sbc); 
[~, ind_fpe]           =    min(fpe);
p_sbc                  =    pMin + ind_sbc - 1;
p_fpe                  =    pMin + ind_fpe - 1

if (strcmp(arfit_opt, 'sbc')) 
    p  = p_sbc;
elseif (strcmp(arfit_opt, 'fpe')) 
    p  = p_fpe;
end



%%    
gamma2  = DTF(ts, lFreq, hFreq, p, Fs);
shuff_gamma2  = DTFsigvalues(ts, lFreq, hFreq, p, Fs, shuffTimes, sigLevel, []);
gamma2  = DTFsigtest(gamma2, shuff_gamma2);
    
  

%%
if( isempty(varargin) )
    rel_err = -1;
else
%     ref_dtf  = varargin{:};
%     for i = 1 : hFreq-lFreq+1
%         ref_total_outflow  = sum(sum(ref_dtf(:, :, i)));
%         gmm_total_outflow  = sum(sum(gamma2(:, :, i)));
%         if (~ref_total_outflow)
%             ref_total_outflow  = 1; 
%         end
%         if (~gmm_total_outflow) 
%             gmm_total_outflow  = 1; 
%         end
%         
%         nlz_ref  = ref_dtf(:, :, i) / ref_total_outflow;
%         nlz_gmm  = gamma2(:, :, i)  / gmm_total_outflow;
%         ref_pow  = norm(nlz_ref, 'fro');
%         if (~ref_pow) 
%             ref_pow = 1; 
%         end
%         rel_err(i)   = norm(nlz_gmm - nlz_ref, 'fro') / ref_pow;
% 
%     end
%     rel_err  = mean(rel_err)
    
    ref_gamma2  = varargin{:};
    ref_dtf  = mean(ref_gamma2, 3);
    ave_gmm  = mean(gamma2, 3);
    ref_total_outflow  = sum(ref_dtf(:));
    gmm_total_outflow  = sum(ave_gmm(:));
    if (~ref_total_outflow)
        ref_total_outflow  = 1; 
    end
    if (~gmm_total_outflow) 
        gmm_total_outflow  = 1; 
    end

    nlz_ref  = ref_dtf / ref_total_outflow;
    nlz_gmm  = ave_gmm / gmm_total_outflow;
    ref_pow  = norm(nlz_ref, 'fro');
    if (~ref_pow) 
        ref_pow = 1;
    end
    rel_err  = norm(nlz_gmm - nlz_ref, 'fro');

end



%%
if( strcmp(dsply_opt, 'Display') )
    h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
    for i = 1 : nnd
        for j = 1 : nnd
            subplot(nnd, nnd, (i-1) * nnd + j)
            plot(lFreq:hFreq, squeeze(gamma2(i, j, :)))
            axis([lFreq, hFreq, 0, 2])
            xlabel('f'), ylabel('\gamma^2')
            title([num2str(j), ' \rightarrow ', num2str(i)])
        end
    end
    
    h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
    image(squeeze(mean(gamma2, 3)), 'CDataMapping','scaled')
    colormap hot
    colorbar
end



end