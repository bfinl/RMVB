function [moment, p_max]  = setup_activity(frq, nnd, nts, con_typ, nlz_opt, dsply_opt, varargin)

if( ~isempty(varargin))
    Spike  = varargin{:};
end

switch con_typ  

        
    case 'Regular1-v1'
        snr  = inf; 
        spike_dur  = 0.080; % second
        spike_pos  = 0.060; % second
        sigma2 = (spike_dur / 4) ^ 2;
        spike_amp  = 10;
        p_max = 0;
        tp = -p_max/frq : 1/frq : (nts-1)/frq;
        moment = zeros(nnd, nts + p_max);
        moment(1, :)  = spike_amp*(exp(-((tp - spike_pos).^2) / sigma2) + ...
                     -0.20*exp(-((tp - spike_pos-spike_dur/3).^2) / sigma2) + ...
                     +0.05*exp(-((tp - spike_pos-spike_dur/3*2).^2) / sigma2));
        moment(1, :)  = moment(1, :) + 10^(-snr/20) * norms(moment(1, :)) / sqrt(size(moment, 2)) * randn(1, nts+p_max);
        moment  = moment(:, p_max+1:end); 
        
        
    case 'Regular1-v2'
        t = 0 : 1/frq : (nts-1)/frq;
        Omega  = 0.09;
        t1  = 0.2;
        t0  = 0.3;
        f  = 13;
        moment(1, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(1, :) = moment(1, :) / max(abs(moment(1, :)));

    
    case 'Regular2-v1'
        snr  = inf; 
        spike_dur  = 0.080; % second
        spike_pos  = 0.060; % second
        sigma2 = (spike_dur / 4) ^ 2;
        spike_amp  = 10;
        p_max = ceil(0.040 * frq);
        tp = -p_max/frq : 1/frq : (nts-1)/frq;
        moment = zeros(nnd, nts + p_max);
        moment(1, :)  = spike_amp*(exp(-((tp - spike_pos).^2) / sigma2) + ...
                     -0.20*exp(-((tp - spike_pos-spike_dur/3).^2) / sigma2) + ...
                     +0.05*exp(-((tp - spike_pos-spike_dur/3*2).^2) / sigma2));
        moment(1, :)  = moment(1, :) + 10^(-snr/20) * norms(moment(1, :)) / sqrt(size(moment, 2)) * randn(1, nts+p_max);
        wgn_amp  = 0.01 * norms(moment(1, :)) / sqrt(size(moment, 2));
        
        for i = p_max + 1 : nts + p_max
            moment(2, i)  = 0.5 * moment(1, i-p_max) + 0.8 * moment(2, i-1) + -0.3 * moment(2, i-2) + wgn_amp * randn(1); 
        end
        moment  = moment(:, p_max+1:end);  
        
        
    case 'Regular3-v1'
        snr  = inf; 
        spike_dur  = 0.080; % second
        spike_pos  = 0.060; % second
        sigma2 = (spike_dur / 4) ^ 2;
        spike_amp  = 10;
        p1  = ceil(0.040 * frq);
        p_max = ceil(0.080 * frq);
        tp = -p_max/frq : 1/frq : (nts-1)/frq;
        moment = zeros(nnd, nts + p_max);
        moment(1, :)  = spike_amp*( ...
                     +1.00*exp(-((tp - spike_pos).^2) / sigma2) + ...
                     -0.20*exp(-((tp - spike_pos-spike_dur/3).^2) / sigma2) + ...
                     +0.05*exp(-((tp - spike_pos-spike_dur/3*2).^2) / sigma2) + ...
                     +1.00*exp(-((tp - 0.25).^2) / sigma2) + ...
                     -0.50*exp(-((tp - 0.30).^2) / sigma2 / 2) + ...
                     +0.40*exp(-((tp - 0.35).^2) / sigma2 / 2) + ...
                     -0.30*exp(-((tp - 0.45).^2) / sigma2 / 4) + ...
                     +0.20*exp(-((tp - 0.50).^2) / sigma2 / 4));
                 
        moment(1, :)  = moment(1, :) + 10^(-snr/20) * norms(moment(1, :)) / sqrt(size(moment, 2)) * randn(1, nts+p_max);
        wgn_amp  = 0.2 * norms(moment(1, :)) / sqrt(size(moment, 2));
        
        for i = p_max + 1 : nts + p_max
            moment(2, i)  = 0.35 * moment(1, i-p1) + 0.25 * moment(2, i-3) +  +0.1 * moment(2, i-5) + wgn_amp * randn(1); 
            moment(3, i)  = 0.6 * moment(1, i-p_max) + 0.1 * moment(3, i-5) + -0.0 * moment(3, i-2) + wgn_amp * randn(1);             
        end
        moment  = moment(:, p_max+1:end); 
        
        
    case 'Regular3-v2'
        snr  = 20; 
        spike_dur  = 0.080; % second
        spike_pos  = 0.060; % second
        sigma2 = (spike_dur / 4) ^ 2;
        spike_amp  = 10;
        p1  = ceil(0.040 * frq);
        p_max = ceil(0.080 * frq);
        tp = -p_max/frq : 1/frq : (nts-1)/frq;
        moment = zeros(nnd, nts + p_max);
        moment(1, :)  = spike_amp*( ...
                     +1.00*exp(-((tp - spike_pos).^2) / sigma2) + ...
                     -0.20*exp(-((tp - spike_pos-spike_dur/3).^2) / sigma2) + ...
                     +0.05*exp(-((tp - spike_pos-spike_dur/3*2).^2) / sigma2) + ...
                     +1.00*exp(-((tp - 0.25).^2) / sigma2) + ...
                     -0.50*exp(-((tp - 0.30).^2) / sigma2 / 2) + ...
                     +0.40*exp(-((tp - 0.35).^2) / sigma2 / 2) + ...
                     -0.30*exp(-((tp - 0.45).^2) / sigma2 / 4) + ...
                     +0.20*exp(-((tp - 0.50).^2) / sigma2 / 4));
                 
        wgn_amp  = 10^(-snr/20) * norms(moment(1, :)) / sqrt(size(moment, 2));         
        moment(1, :)  = moment(1, :) + wgn_amp * randn(1, nts+p_max);
        
        for i = p_max + 1 : nts + p_max
            moment(2, i)  = 0.35 * moment(1, i-p1) + 0.25 * moment(2, i-3) +  +0.1 * moment(2, i-5) + wgn_amp * randn(1); 
            moment(3, i)  = 0.6 * moment(1, i-p_max) + 0.1 * moment(3, i-5) + -0.0 * moment(3, i-2) + wgn_amp * randn(1);             
        end
        moment  = moment(:, p_max+1:end); 
        
        
    case 'Regular3-v3'
        Spk_ts = [Spike(59:170,66);zeros(88,1)];
        Spk_ts = Spk_ts/max(Spk_ts);
        dur_len = size(Spk_ts,1);
        tp = [1:dur_len]/frq;
        x1 = Spk_ts + 0.05*randn(dur_len,1);
        x_1_2 = 0.05*randn(dur_len,1);
        x_1_3 = 0.05*randn(dur_len,1);
        for i_1 = 26:dur_len
            x_1_2 (i_1) = 0.35*x1(i_1-25) + 0.1*x_1_2(i_1-5) + + 0.25*x_1_2(i_1-3) + 0.05*randn;
            if i_1 > 20
                x_1_3 (i_1) = 0.9*x_1_2(i_1-20) + (0.1)*x_1_3(i_1-5)+ 0.05*randn;
            elseif i_1 > 30
            x_1_3 (i_1) = 0.5*x1(i_1-30) + 0.9*x_1_2(i_1-20) + (0.1)*x_1_3(i_1-5)+ 0.05*randn;
            end
        end
        moment(1, :) = x1;
        moment(2, :) = x_1_2;
        moment(3, :) = x_1_3;          
        
        
    case 'Regular3-v4'
        t = 0 : 1/frq : (nts-1)/frq;
        Omega  = 0.09;
        t1  = 0.15;
        t0  = 0.05;
        f  = 9;
        moment(1, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(1, :) = moment(1, :) / max(abs(moment(1, :)));
        
        Omega  = 0.09;
        t1  = 0.2;
        t0  = 0.3;
        f  = 13;
        moment(2, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(2, :) = moment(2, :) / max(abs(moment(2, :)));
        
        Omega  = 0.05;
        t1  = 0.3;
        t0  = 0.1;
        f  = 11;
        moment(3, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(3, :) = moment(3, :) / max(abs(moment(3, :)));
        
        
    case 'Regular3-v5'
        t = 0 : 1/frq : (nts-1)/frq;
        Omega  = 0.09;
        t1  = 0.15;
        t0  = 0.05;
        f  = 9;
        moment(1, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(1, :) = moment(1, :) / max(abs(moment(1, :)));
        
        Omega  = 0.09;
        t1  = 0.2;
        t0  = 0.3;
        f  = 13;
        moment(2, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(2, :) = moment(2, :) / max(abs(moment(2, :)));
        
        Omega  = 0.05;
        t1  = 0.3;
        t0  = 0.1;
        f  = 11;
        moment(3, :) = exp(-((t-t1)./Omega).^2) .* sin(2*pi*f*(t-t0));
        moment(3, :) = moment(3, :) / max(abs(moment(3, :)));
        
        moment(2, :) = 0.5 * moment(1, :) + 0.5 * moment(2, :);
        moment(3, :) = 0.5 * moment(2, :) + 0.5 * moment(3, :);
        
        
    case 'Extent'
         
        
    otherwise
        display('*************** Weird Request ***************')
        
        
end



%%

if (strcmp(nlz_opt, 'Normalize')) 
    moment  = (moment - repmat(mean(moment, 2), [1, nts])) ./ repmat(std(moment, 1, 2), [1, nts]);
end

        
%%

if( strcmp(dsply_opt, 'Display') )
    
    for i = 1 : nnd lgnd{i} = ['Source ' num2str(i)]; end
    tp = 0:1/frq:(nts-1)/frq; tp = repmat(tp.', 1, nnd); tp = tp * 1000;
    h = figure;  ApplyProperties(h, 'Customized-v1'); 
    plot(tp, (randn(1, length(tp))).'), legend(lgnd, 'Location', 'northwest'), xlabel('Time (ms)'), ylabel('Moment (\muAmm)')
    
    h = figure;  ApplyProperties(h, 'Customized-v1'); 
    for i = 1 : nnd 
        lgnd = ['Source ' num2str(i)]; 
        fp = linspace(0, frq*(nts-1)/nts, nts); 
        fp = [fp(ceil(nts/2)+1:end)-frq fp(1:ceil(nts/2))];
        subplot(nnd, 1, i), plot(fp, abs(fftshift(fft(moment(i, :))))) 
        legend(lgnd, 'Location', 'northwest'), xlabel('Frequency (Hz)'), ylabel('Magnitude')
    end
    
%     h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
%     for i = 1 : nnd
%         for j = 1 : nnd
%             subplot(nnd, nnd, (i-1) * nnd + j)
%             mu  = abs(moment(i, :) * moment(j, :).' / norm(moment(i, :)) / norm(moment(j, :)));
%             plot([1 2], mu * [1 1])
%             ylim([0 1+1e-5])
%             ylabel('CF')
%             title([num2str(j), ' \rightarrow ', num2str(i)])
%         end
%     end
    
%     h = figure;  h = ApplyProperties(h, 'Customized-v1'); 
%     for i = 1 : nnd
%         for j = 1 : nnd
%             C(i, j)  = abs(moment(i, :) * moment(j, :).' / norm(moment(i, :)) / norm(moment(j, :)));
%         end
%     end
%     colormap hot, imagesc(C, [0, 1]), colorbar, title('Correlation Factor')
%     h = gca; h.XTick = 1:nnd; h.YTick = 1:nnd;
    
end

end


