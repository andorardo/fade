function res = gaussianBlur(s,sigma)
    res = s;
    if sigma>0
        n = ceil(4*sigma);
        sigma2 = sigma*sigma;
        w = exp(-0.5*(-n:n).^2/sigma2) / sqrt(2*pi*sigma2);
        w = w/sum(w);
        res = filter(w, 1, [s(:); zeros(2*n,1)]);
        res = res(n+1:end-n);
        assert( length(res) == length(s) );
        res = reshape(res,size(s));
    end
end