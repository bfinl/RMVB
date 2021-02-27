function [Phi, pow_inf, src_activity] = add_noise(Phi, type, snr, bll, K, src_activity, lambda, acr)

[noe, nts]  = size(Phi);
Wnoise  = zeros(noe, nts+bll);
Cnoise  = zeros(noe, nts+bll);
awgn_var  = inf;
corr_var  = inf;
awgn_amp  = inf;
corr_amp  = inf;


%%
sig_var  = var(Phi(:), 1);
src_noi_var  = 0;
if( isempty(acr) ), acr  = 1:size(Phi, 2); end
sig_max  = max(norms(Phi(:, acr))) / sqrt(size(Phi, 1));
sig_amp  = norm(Phi(:, acr), 'fro') / sqrt(numel(Phi(:, acr)));


%%
switch type
    
    case 'AWGN'
        
        Wnoise = randn(noe, nts+bll) * sig_amp * (10 ^ (-snr / 20));
        awgn_var  = var(Wnoise(:), 1);
        awgn_amp  = norm(Wnoise, 'fro') / sqrt(numel(Wnoise));
        
                
    case 'ACGN'
        
        % Accounts for both awgn and the interference comming form activities in the brain 
        if( isempty(src_activity) ) 
            src_activity  = randn(size(K, 2), nts+bll); 
            Cnoise  = K * src_activity; 
            const   = 1 / norm(Cnoise, 'fro') * sqrt(numel(Cnoise)) * sig_amp * (10 ^ (-snr / 20)) * sqrt(1-lambda); 
            src_activity  = src_activity * const;
            Cnoise  = Cnoise * const;
            src_noi_var  = const^2; 
            Wnoise  = randn(noe, nts+bll) * sig_amp * (10 ^ (-snr / 20)) * sqrt(lambda);
        else
            Cnoise  = K * src_activity; 
            wnoise_amp  = sqrt((sig_amp * (10 ^ (-snr / 20)))^2  -  (norm(Cnoise, 'fro') / sqrt(numel(Cnoise)))^2);
            if( ~isreal(wnoise_amp) ) 
                wnoise_amp = 0;
                disp('Warning: SNR may not be achieved!!!')
            end
            Wnoise  = randn(noe, nts+bll) * wnoise_amp;
        end
        awgn_var  = var(Wnoise(:), 1);
        corr_var  = var(Cnoise(:), 1);
        awgn_amp  = norm(Wnoise, 'fro') / sqrt(numel(Wnoise));
        corr_amp  = norm(Cnoise, 'fro') / sqrt(numel(Cnoise));

        
    otherwise
        display('*************** Weird Request ***************')


end


%%
noise  = Wnoise + Cnoise;
Phi  = [zeros(noe, bll) Phi];
Phi  = Phi + noise;

pow_inf(1, 1:3)  = [sig_var, awgn_var, corr_var];
pow_inf(2, 1:3)  = [sig_amp, awgn_amp, corr_amp] .^ 2;
pow_inf(3, 1:3)  = [0, 0, src_noi_var];

disp_opt  = 'No Display';
if( strcmp(disp_opt, 'Display') )

    snr_awgn  = 20 * log10(sig_amp / norm(Wnoise, 'fro') * sqrt(numel(Wnoise)));
    snr_corr  = 20 * log10(sig_amp / norm(Cnoise, 'fro') * sqrt(numel(Cnoise)));
    snr_cbnd  = 20 * log10(sig_amp / norm(noise , 'fro') * sqrt(numel(noise)));
    disp('****************************************************')
    fprintf('snr-awgn: %f\nsnr-corr: %f\nsnr-cbnd: %f\n', snr_awgn, snr_corr, snr_cbnd)
    disp('****************************************************')
    
    psnr_awgn  = 20 * log10(sig_max / norm(Wnoise, 'fro') * sqrt(numel(Wnoise)));
    psnr_corr  = 20 * log10(sig_max / norm(Cnoise, 'fro') * sqrt(numel(Cnoise)));
    psnr_cbnd  = 20 * log10(sig_max / norm(noise , 'fro') * sqrt(numel(noise)));
    disp('****************************************************')
    fprintf('psnr-awgn: %f\npsnr-corr: %f\npsnr-cbnd: %f\n', psnr_awgn, psnr_corr, psnr_cbnd)
    disp('****************************************************')
    
    figure, plot(Phi.')
    
end

end

