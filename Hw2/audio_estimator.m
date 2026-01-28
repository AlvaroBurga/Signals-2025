function [R_right,R_left] = audio_estimator(maxLagTime, M, Fs, source, mu, mic_left, mic_right)
    R_right = zeros(M,1);
    R_left = zeros(M,1);
    
    % Parameters
    K1 = floor((maxLagTime*4) * Fs) + 1;
    K2 = K1;
    
    a_left = zeros(K1 + K2+ 1,1);
    a_right = zeros(K1 + K2+ 1,1);
    % 1. Compute tau_i only every K samples
    for i = K1*2:M-2*K2
        %tau = tau_estimator(t(i), mic_left, mic_right, Fs, maxLagTime);
        %display_index = floor(tau * Fs);
        display_index = 0;
        rx_right = mic_right(i - K1 + display_index: i + K2 + display_index, 1);
        rx_left = mic_left(i - K1: i + K2);
        
        R_right(i) = a_right'*rx_right;
        R_left(i) = a_left' * rx_left;
    
        eps_right= source(i) - R_right(i);
        eps_left = source(i) - R_left(i);
    
        trR_r = rx_right'*rx_right;
        trR_l = rx_left'*rx_left;
    
        a_right=a_right+mu*eps_right*rx_right/trR_r;
        a_left=a_left+mu*eps_left*rx_left/trR_l;
    end

end