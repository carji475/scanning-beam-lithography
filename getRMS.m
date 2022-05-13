function rms_error = getRMS(litoset, result)
fft_exposure(result.Nxm, result.Nym, result.support_x, result.support_y, ...
    result.xres, result.yres, result.Hsqn);
XX = fft_exposure(result.wwTot); % Dosage
fx = 1./( 1+ exp(-litoset.a*(XX-litoset.tr)) ); % feature
rms_error = rms(litoset.Zm(:)-fx(:));
end

