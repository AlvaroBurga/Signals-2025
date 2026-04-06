function [mic,ToD,IoD] = rx_generator(SNR, Fs, x, x_source, y_source,vp, s, i1, i2)

    M = i2 - i1 + 1; %Number of samples
    SNR_linear = 10^(SNR/10); 
    var_w = 1/SNR_linear;
    w = sqrt(var_w)*randn(M,1); %generated noise
    r = sqrt((x-x_source).^2 + y_source.^2); %Distance of the link
    ToD = r/vp; %Time of delay
    IoD = floor(ToD*Fs) + 1; %Index of delay
    idx  = (i1:i2) + IoD(:)'; 
    s_delayed  = s(idx , :); %Delaying the signal
    mic = attenuation(r).*(s_delayed(:,1) + s_delayed(:,2)) + w; %Recieved signal in the mic
end