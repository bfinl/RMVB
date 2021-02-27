function [LFD, gammamat] = add_uncertainty(LFD, snr, lfd_type)

if(strcmp(lfd_type, 'rot')), d = 3; elseif(strcmp(lfd_type, 'fxd')), d = 1; end

m = size(LFD, 1);
gammamat  = zeros(m, m, size(LFD, 2)/d);
for i = 1 : d : size(LFD, 2)
    sigma = norm(LFD(:, i:i+d-1), 'fro') / sqrt(m * d) / (10^(snr/20)) / sqrt(d);
    LFD(:, i:i+d-1) = LFD(:, i:i+d-1) + sigma * randn(m, d);
    gammamat(:, :, ceil(i/d)) = sqrt(chi2inv(0.975, m)) * sigma * eye(m);
end

end