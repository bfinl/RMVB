function FilteredPhi = my_filter(Phi, frq, cutofffrq, filt_type, bll, disp_opt)

nts  = size(Phi, 2);
tps  = linspace(-bll/frq, (nts-bll-1)/frq, nts);
fps  = linspace(0, frq*(nts-1)/nts, nts); % frequency points
sfp  = [fps(ceil(nts/2)+1:end)-frq fps(1:ceil(nts/2))]; % shifted frequency points   
FilteredPhi = zeros(size(Phi));
lcutoff  = cutofffrq(1);
hcutoff  = cutofffrq(2);

switch filt_type
    
    case 'FiltFilt'
        b = fir1(ceil(0.1*nts), [lcutoff hcutoff]/(frq/2));
%         b = fir1(floor(size(Phi, 2)/3), [lcutoff hcutoff]/(frq/2)]);
        for i = 1 : size(Phi, 1)
            FilteredPhi(i, :)  = filtfilt(b, 1, Phi(i, :));
        end
        if(strcmp(disp_opt, 'Display'))
            h = figure; ApplyProperties(h, 'Customized-v4'); freqz(b, 1);
        end
        
    case 'FFTiFFT'
        for i = 1 : size(Phi, 1) 
            filtermap  = (fps>=lcutoff & fps<=hcutoff) | (fps<=frq-lcutoff & fps>=frq-hcutoff);
            FilteredPhi(i, :) = real(ifft(fft(Phi(i, :)) .* filtermap));
        end
        
    case 'EEGLab'
        EEG.trials = 1;
        EEG.event  = [];
        EEG.srate  = frq;
        EEG.data   = Phi;
        EEG.pnts   = nts;
        [EEG, ~, b]  = pop_eegfiltnew(EEG, lcutoff, hcutoff, [], 0, [], strcmp(disp_opt, 'Display'), []);
        fprintf('Filter Order: %f\n', length(b)-1)
        FilteredPhi  = EEG.data;
        
    case 'NULL' 
        FilteredPhi  = Phi;
 
end


%%
if(strcmp(disp_opt, 'Display'))

    h = figure; ApplyProperties(h, 'Customized-v4');
    plot(tps, Phi.'), title('Noisy Data'), xlabel('Time (ms)'), ylabel('Amplitude (v)')
    h = figure; ApplyProperties(h, 'Customized-v4'); 
    plot(tps, FilteredPhi.'), title('Filtered Data'), xlabel('Time (ms)'), ylabel('Amplitude (v)')

    
    h = figure; ApplyProperties(h, 'Customized-v4');
    plot(sfp, abs(fftshift(fft(Phi, [], 2), 2)).'), title('Noisy Data'), xlabel('Frequency (Hz)'), ylabel('Magnitude (v)')
    h = figure; ApplyProperties(h, 'Customized-v4');
    plot(sfp, abs(fftshift(fft(FilteredPhi, [], 2), 2)).'), title('Filtered Data'), xlabel('Frequency (Hz)'), ylabel('Magnitude (v)')
    
end

end