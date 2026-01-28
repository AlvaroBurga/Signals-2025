%REVIEW APENDIX 20B FOR THE THEORY
function delta_tau = tau_estimator(t0, mic_left, mic_right, Fs, maxLagTime, resolution)
    T = maxLagTime*4; %Use the maxlag arround t0
    K1 = max(1, floor((t0-T/2) * Fs) + 1);
    K2 = min(length(mic_right),floor((t0+T/2) * Fs) + 1);

    x = mic_left(K1:K2);
    y = mic_right(K1:K2);

    max_lag = floor(maxLagTime * Fs);

    [r, lags] = xcorr(x, y, max_lag, 'coeff');

    [~, idx0] = max(r);
    lag0 = lags(idx0);
    
    if resolution > 0
        N = length(r);
    
        % higher tau resolution
        tau_min = lag0 - 1;     % search ±1 sample
        tau_max = lag0 + 1;
        tau_vec = linspace(tau_min, tau_max, resolution);
    
        % Function h(t)
        h = @(t) sin(pi*t) ./ ( N * sin(pi*t/N) + eps );   % +eps avoids division by zero
    
        % ========== Evaluate interpolated correlation ===================
        r_interp = zeros(size(tau_vec));
        
        for i = 1:length(tau_vec)
            tau = tau_vec(i);
            a = h( tau - lags );
            r_interp(i) = sum( r .* h( tau - lags )' );
        end
    
        % === Find the interpolated maximum and convert to time ===
        [~, idx_up] = max(r_interp);
        delta_index = tau_vec(idx_up); %Delta index is fractional. According with the resolution
        delta_tau = delta_index / Fs;
    else 
        delta_tau = lag0 / Fs; % Use the lag0 value for delta_tau if interpolation is not performed
    end
    
end