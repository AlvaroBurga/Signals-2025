function [mic,ToD,IoD] = rx_generator_no_noise( Fs, x, x_source, y_source,vp, s, i1, i2)

r = sqrt((x-x_source).^2 + y_source.^2); %Distance of the link
ToD = r/vp; %Time of delay
IoD = floor(ToD*Fs) + 1; %Index of delay
idx  = (i1:i2) + IoD(:)'; 
s_delayed  = s(idx , :); %Delaying the signal
mic = attenuation(r).*(s_delayed(:,1) + s_delayed(:,2)); %Recieved signal in the mic
end