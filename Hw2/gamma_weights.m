function w = gamma_weights(N, sigma)
    % N     : vector lenght
    % sigma : dacaying factor 
    
    center = (N + 1) / 2;
    x = 1:N;
    
    w = exp(-0.5 * ((x - center) / sigma).^2);
    w = w / sum(w);  
end