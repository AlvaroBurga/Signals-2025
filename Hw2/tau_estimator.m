%REVIEW APENDIX 20B FOR THE THEORY
function delta_tau = tau_estimator(t0, mic_left, mic_right, Fs, win_size, resolution)
    T = win_size / Fs; % Calculate the time window in seconds
    K1 = max(1, floor((t0-T/2) * Fs) + 1); %Start point of the window
    K2 = min(length(mic_right),floor((t0+T/2) * Fs) + 1); %end point of the window

    left = mic_left(K1:K2); %left data used
    right = mic_right(K1:K2); %right data used

    [r, lags] = xcorr(left, right, win_size, 'coeff'); %Calculating the correlation in the window

    [~, idx0] = max(r.^2); %Obtaining the maximun correlation
    lag0 = lags(idx0); %It's index
    
    if resolution > 0
        N = length(r);
    
        % higher tau resolution
        tau_min = lag0 - 1;     % search ±1 sample
        tau_max = lag0 + 1;
        tau_vec = linspace(tau_min, tau_max, resolution);
    
        % Function h(t)
        h = @(t) sin(pi*t) ./ ( N * sin(pi*t/N) + eps );   % +eps avoids division by zero
    
        % Evaluate interpolated correlation
        r_interp = zeros(size(tau_vec));
        
        for i = 1:length(tau_vec)
            tau = tau_vec(i);
            r_interp(i) = sum( r .* h( tau - lags )' );
        end
    
        % Find the interpolated maximum and convert to time
        [~, idx_up] = max(r_interp);
        delta_index = tau_vec(idx_up); %Delta index is fractional. According with the resolution
        delta_tau = delta_index / Fs;
    else 
        delta_tau = lag0 / Fs; % Use the lag0 value for delta_tau if interpolation is not performed
    end
    
end