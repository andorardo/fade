function res = gaussianDetrend(data, sampling, fcut, lowPass)
    if nargin<4, lowPass = false; end
    n = length(data);
    fftdata = reshape(fft(data),[],1);
    T = n/sampling;
    if length(fcut)==1
        edlen = floor(min(15*T*fcut, n/2));
        if lowPass, edlen = n/2; end
        w = 1.0 - exp(-0.5*((0:edlen)'/(T*fcut)).^2);
        if lowPass, w = 1.0-w; end
        fftdata(1:edlen+1) = fftdata(1:edlen+1) .* w;
        fftdata(end-edlen+1:end) = conj(fftdata(edlen+1:-1:2));
    else
        error('Not yet handling band-pass filtering');
    end
    res = reshape(ifft(fftdata),size(data));
end