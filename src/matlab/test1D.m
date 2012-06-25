function m = test1D(xdim,fdim)
    % Test deconvolution inference
    if nargin<1
        xdim = 100;
    end
    if nargin<2
        %fdim = ceil(sqrt(xdim-1));
        fdim = ceil(10);
        fdim = fdim+mod(fdim,2)-1;
    end
    assert(mod(fdim,2)==1);

    sorig = zeros(xdim,1);
    p = randperm(xdim);
    nSpikes = ceil(0.1*xdim);
    sorig(p(1:nSpikes)) = rand(nSpikes,1);
    f = zeros(2*fdim,1);
    f(1:fdim) = normpdf(6*((0:fdim-1)/(fdim-1)-0.5),0,1)/normpdf(0,0,1);
    x = conv2(f,sorig);
    x = x(1:xdim);
    sigma2 = (0.03)^2;
    x = x+sqrt(sigma2)*randn(size(x));

    figure(1); clf;
    subplot(5,4,1);
    pfun(sorig);
    title('sorig');
    subplot(5,4,2);
    pfun(f);
    title('f');
    subplot(5,4,3);
    pfun(x);
    title('x');
    drawnow;

    m = MEMDeconvolution(x,f);
    m.sigma2 = sigma2;
    m.alpha = 1e3;
    m.stopMethod = 1; % Autonoise
    m.MaximizeEvidenceAlphaOnly;
    subplot(5,4,4);
    s = reshape(m.s,xdim,1);
    pfun2(s);
    title(['s, alpha=' num2str(m.alpha)]);
    drawnow;

    alphas = logspace(3,-4,8);
    for k=1:length(alphas)
        m.alpha = alphas(k);
        m.sConverge;
        subplot(5,4,4+2*k-1);
        s = reshape(m.s,xdim,1);
        pfun2(s);
        title(['alpha=' num2str(m.alpha)]);
        drawnow;
        subplot(5,4,4+2*k);
        pfun(s-sorig);
        title(['norm=' num2str(norm(s(:)-sorig(:)))]);
        drawnow;
    end

    function pfun2(x)
        stairs(1:length(x(:)),x(:),'b-');
        hold on; stairs(1:length(sorig(:)),sorig(:),'r--');
    end
    function pfun(x)
        stairs(1:length(x(:)),x(:),'b-');
    end

end
