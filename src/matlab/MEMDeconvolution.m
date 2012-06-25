classdef MEMDeconvolution < handle
%Deconvolution using quantified maximum entropy prior
%   
%   Accompanies the following article:
%
%   Fast, automated implementation of temporally precise blind deconvolution of multiphasic excitatory postsynaptic currents
%
%   Daniel Andor-Ardó (1)(2), Erica C. Keen (1), A. J. Hudspeth (1), Marcelo O. Magnasco (2)
%
%   (1) Howard Hughes Medical Institute and Laboratory of Sensory Neuroscience, and
%   (2) Laboratory of Mathematical Physics,
%   The Rockefeller University, New York, USA
%
%   PLoS ONE 2012
%   http://dx.doi.org/10.1371/journal.pone.0038198
%
%   Author: Daniel Andor-Ardó, daniel.andor@gmail.com

    
    properties (Access=public)
        x; xr; xc; % rows and columns
        s;
        f; fr; fc;
        regMethod; % 1==MaxEnt[flux]
        stopMethod = 2; % 1==Classic; 2==ClassicAutonoise
        threshold = 5;
        smoothing = 0;
        alpha   = 100. ;
        sigma2  = 1. ;
        
        Gfactor = 2. ;
        smaxiter = 100;
        
        % for line search
        BStol = 0.01;
        BSsteptol = 1e-5;
        BSmaxrec = 6;
        LMlastx = 0.01;
        LMgtyp = 1.0;
        LMmaxits = 25*2*5/2;
        LMstep = 1.0;
        LMpath = {};
    end
    
    properties (Access=private)
        % Cache, set by setters
        FTx;
        FTFx;
        
        % Cache, handled explicitly
        storedDsdU = [];
        storedG = [];
        r = [];
    end
    
    properties (Dependent=true)
        beta;
        sparseS;
        reconstruction;
        residual;
    end
    
    methods

        function obj = set.x(obj,x_) %#ok<*MCHV2>
            obj.x = x_(:);
            [obj.xr,obj.xc] = size(x_);
            if ~isempty(obj.f)
                obj.init;
            end
            % Now possible set other stuff, see SetX[]
        end
        
        function obj = set.f(obj,f_)
            obj.f = f_(:);
            [obj.fr,obj.fc] = size(f_);
            if ~isempty(obj.x)
                obj.init;
            end
        end
        
        function b = get.beta(obj)
            b = 1./obj.sigma2;
        end
        
        function x = get.reconstruction(self)
            x = self.ApplyF(self.sparseS);
        end
        
        function r = get.residual(self)
            r = self.x - self.reconstruction;
        end
        
        function sparseS = get.sparseS(self)
            a = self.amplitudes;
            t = self.times;
            sparseS = zeros(size(self.s));
            sparseS(t) = a;
        end
        
    end

    methods (Access=public)
        
        function obj = MEMDeconvolution(x_,f_,regMethod_,threshold_,smoothing_)
            if nargin<3 || isempty(regMethod_)
                regMethod_ = [1 mean(x_(:))/sum(f_(:))]; % MaxEnt[flux] -- better calc flux
            end
            obj.regMethod = regMethod_;

            if nargin>3 && ~isempty(threshold_)
                obj.threshold = threshold_;
            end
            if nargin>4 && ~isempty(smoothing_)
                obj.smoothing = smoothing_;
            end
            
            obj.f = f_;
            obj.x = x_;
            obj.resetS;
        end
        
        function sConverge(self)
            iter = 1;
            oldXGratio = 1.;
            Alpha = 1; %NRDampingAlpha
            self.storedDsdU = [] ;
            self.storedG = [] ;
            X = self.CrossEntropy(self.s);
            G = self.ApproximateGalpha(self.s, 1); 
            while X >= G/(10*self.Gfactor) ...
                    && (oldXGratio/(X/G) - 1.)^2 >= 0.01 ...
                    && iter < self.smaxiter
                iter = iter + 1;
                oldXGratio = X/G;
                [Deltas, g] = self.GetDsdU(self.s);
                self.s = self.s + Alpha * Deltas;
                self.storedDsdU = [] ;
                self.storedG = [] ;
                switch self.regMethod(1)
                    case 1 %'MaxEnt[_] | MaxEnt2 | PositivityPrior[_] | Shannon'
                        self.s = max(self.s,sqrt(eps));
                    otherwise %'WienerFilter | EmptinessPrior[_, _], '
                        % nada
                end
                X = self.CrossEntropy(self.s);
                G = self.ApproximateGalpha(self.s, 1); 
            end
        end
        
        function MaximizeEvidenceAlphaOnly(self,varargin)
            bisect = false;
%             x = log(self.alpha); % params to optimize
%             MEpath = {};
            self.sConverge;
            grad1 = -self.gradLogE1();
            dir = -grad1;
            self.LMReset();
%             LMgtyp = sqrt((grad'*grad)/length(self.x));
            x = self.LineMinimize1(dir,bisect); % Use bisection for more precision
            self.alpha = exp(x);
            self.sConverge;
%             newgrad = -self.gradLogE1();
%             MEpath{end+1} = [dir'*grad,X  XXX LMpath, {{LMlastx/LMgtyp, dir.newgrad}}]];
        end
        
        function LMReset(self)
            self.BStol = 0.001;
            self.BSsteptol = 10^-5;
            self.BSmaxrec = 6;
            self.LMlastx = 0.01;
            self.LMgtyp = 1.0;
            self.LMmaxits = 25*2*5/2;
            self.LMstep = 1.0;
            self.LMpath = {};
        end
        
        function alphas = getSortedAlphaLMPath(self)
            alphas = cellfun(@(c) c{3}, self.LMpath);
            [~,IX] = sort(alphas);
            IX = IX(alphas(IX)>=alphas(end));
            alphas = alphas(IX(end:-1:1));
        end
        
        function newparams = LineMinimize1(self,dir,doBisection)
            % In alpha direction (1) *only*
            x = self.LMlastx / self.LMgtyp;
            origparams = log(self.alpha);
	
            function [dirggx,ggx] = prod(x)
                newp = origparams + x * dir;
                self.alpha = exp(newp);
                self.sConverge;
                ggx = -self.gradLogE1;
                self.LMpath{end+1} = {x, dir*ggx, self.alpha}; 
                dirggx = dir'*ggx;
            end
            
            function res = bisection(a, b, fa, fb, depth)
                c = (a + b)/2;
                %if min(abs(fa),abs(fb))/length(self.x) < self.BStol || depth > self.BSmaxrec
                if abs(fa-fb)/max(abs(fa),abs(fb))/length(self.x) < self.BStol || depth > self.BSmaxrec
                    ma=abs(fa);
                    mb=abs(fb);
                    c = ma+mb;
                    ma = ma/c; mb = mb/c;
                    res = mb*a + ma*b; % (* linear interpolation *)
                else
                    [fc,~] = prod(c); %#ok<*PROP>
                    logInfo;
                    if fc*fa < 0
                        res = bisection(a, c, fa, fc, depth+1);
                    else
                        res = bisection(c, b, fc, fb, depth+1);
                    end
                end
            end
            
            function logInfo
                fprintf('its=%d alpha=%f\n',its,self.alpha);
            end
            
            scaleFactor = sqrt(10);
            LMmaxalphascale = 1/scaleFactor;
            
            x = min(-log(LMmaxalphascale)/abs(dir), x);
            [ms, gx] = prod(x);
            its = 1;
            if ms < 0 % (* Go further *)
                while its <= self.LMmaxits
                    logInfo;
                    y = scaleFactor * x;
                    y = x + min(-log(LMmaxalphascale)/abs(dir), y-x);
                    [mt, gy] = prod(y);
                    if mt >= 0, break; end
                    x = y; ms = mt; tmp=gx; gx=gy; gy=tmp;
                    its = its+1;
                end
            elseif ms > 0 %  (* Step back *)
                while its <= self.LMmaxits
                    logInfo;
                    y =  x / scaleFactor;
                    [mt, gy] = prod(y);
                    if mt <= 0, break; end
                    x = y; ms = mt; tmp=gx; gx=gy; gy=tmp;
                    its = its+1;
                end
            else % s = 0, 'hole in one'
                t = 1.0; y = x;
            end
            ms = abs(ms);
            mt = abs(mt);
            mm = ms+mt;
            ms = ms/mm; mt = mt/mm;
            mm = ms*y + mt*x; 

            if doBisection
                if y > x
                    mm = bisection(x, y, gx'*dir, gy'*dir, 1);
                else
                    mm = bisection(y, x, gy'*dir, gx'*dir, 1);
                end
            end
            self.LMlastx = mm*1.0*self.LMgtyp;
            newparams = origparams + mm*dir;
            self.LMstep = norm(mm*dir)/3;
        end
        
        function alphaConverge(self,alphaStart,maxIter)
            % Deprecated in favor of MaximizeEvidenceAlphaOnly
            if nargin<2
                alphaStart = self.alpha;
            end
            if nargin<3
                maxIter = 100;
            end
            iter = 1;
            q=[];
            Lqs = [];
            Lalphas = [];
            tol = 0.0001;
            
            self.alpha = alphaStart;
            fprintf('iter=%d alpha=%f\n',iter,self.alpha);
            Lalphas(end+1) = log2(self.alpha);
            self.sConverge;
            q = self.ApproximateQ(self.s); %XX
            Lqs(end+1) = log2(q); 

            if Lqs(iter)^2 < tol, return; end
            iter = iter+1;
            if q < 1.
                self.alpha = self.alpha/(1.1+rand*0.1);
            else
                self.alpha = self.alpha*(1.1+rand*0.1);
            end
            fprintf('iter=%d alpha=%f\n',iter,self.alpha);
            Lalphas(end+1) = log2(self.alpha);
            self.sConverge;
            q = self.ApproximateQ(self.s); 
            Lqs(end+1) = log2(q);

            da0 = log2(10);
            da = da0;
            while Lqs(iter)^2 > tol && iter < maxIter
                if diff(Lalphas(end-1:end))>0
                    da = da/2;
                else
                    da = min(da*1.1,da0);
                end
                deltaLalpha = -Lqs(iter)*(Lalphas(iter) - Lalphas(iter-1))/(Lqs(iter) - Lqs(iter-1));
                deltaLalpha = min(da, max(-da, deltaLalpha));
                Lalphas(end+1) = log2(self.alpha) + deltaLalpha;
                iter = iter + 1;
                self.alpha = 2^Lalphas(iter);
                self.storedDsdU = [];
                self.storedG = []; 
                fprintf('iter=%d alpha=%f\n',iter,self.alpha);
                self.sConverge; 
                q = self.ApproximateQ(self.s); 
                Lqs(end+1) = log2(q);
            end
        end

        function a = amplitudes(self)
            a = sumSAmplitudes(self.s,self.threshold,self.smoothing);
        end
        
        function t = times(self)
            t = peakSAmplitudes(self.s,self.threshold,self.smoothing);
        end
        
        function f = filterEstimate(self, varargin)
            % Uses x and s to estimate f

            % Optional arguments; ls==LeastSquares, lp==LinearProgramming
            method = 'ls';
            if ~isempty(which('linprog')) % check to see if Optimization Toolbox function 'linprog' is available
                method = 'lp';
            end
            flen = length(self.f);
            peakPosition = 25;
            refinePeak = true;
            varlen = length(varargin);
            i_ = 0;
            while i_<varlen
                i_ = i_ + 1;
                switch lower(varargin{i_})
                    case 'method'
                        i_ = i_ + 1;
                        method = varargin{i_};
                    case 'filterlength'
                        i_ = i_ + 1;
                        flen = varargin{i_};
                    case 'peakposition'
                        i_ = i_ + 1;
                        peakPosition = varargin{i_};
                    case 'refinepeak'
                        i_ = i_ + 1;
                        refinePeak = varargin{i_};
                end
            end
            
            sparseS = self.sparseS;
            
            sSS = fft(sparseS);
            XX = fft(self.x);
            Cax1 = ifft(conj(sSS) .* XX);
            Caa1 = ifft(sSS .* conj(sSS));
            tmat = abs(bsxfun(@minus, (1:flen), (1:flen)')) + 1;
            m = Caa1(tmat);
            b = Cax1(1:flen);

            switch lower(method)
                case 'ls'
                    f = m \ b;
                case 'lp'
                    f = linprog(ones(size(b)),-m,-b);
                otherwise
                    error(['Unknown method ''' method '''']);
            end
            [~,I] = max(f);
            f = f-min(f(1:I));
            assert( size(f,2)==1 );
            delta = peakPosition-I;
            if delta>0
                f = [f(1)*ones(delta,1);f(1:end-delta)];
            end
            if delta<0
                f = [f(abs(delta)+1:end);f(end)*ones(abs(delta),1)];
            end
            f = f-max(0,min(f(1:peakPosition)));
            f = f/max(f);
            
            if refinePeak
                %NX = length(self.x);
                NX = 2*flen;
                iomegas = 2i*pi*[0:NX/2-1, -NX/2:-1]';
                fpadded = f(end)*ones(NX,1);
                fpadded(1:flen) = f;
                FF = fft(fpadded);
                df = real(ifft(iomegas.*FF));
                tmp = fzero(@(t) interp1(1:length(df),df,peakPosition+t,'spline'), 0);
                f = real(ifft(FF.*exp(-iomegas*tmp/NX)));
                f = f(1:flen);
            end
        end

    end
    
    methods (Access=protected)
        
        function init(obj,resetS)
            obj.FTx = obj.ApplyFT(obj.x); %#ok<*MCSUP>
            obj.FTFx = obj.ApplyFTF(obj.x);
        end
        
        function resetS(obj)
            obj.s = ones(size(obj.x))/length(obj.x);
            if length(obj.regMethod)>1
                obj.s = ones(size(obj.x))*obj.regMethod(2);
            end
        end
        
    end
    
    methods (Access=private)
        
        function q = ApproximateQ(self,s)
            switch self.stopMethod
                case 1 % Classic
                    q = self.ApproximateGalpha(s)/(2*self.alpha*self.H(s));
                case 2 % ClassicAutonoise
                    q = self.U(s) * self.ApproximateGalpha(s)/(self.alpha * self.H(s) * length(self.x));
            end
        end
        
        function res = gradLogE1(self)
            % alpha only
            % w.r.t. log(alpha): use df/dlog(a) = a df/da.
    		res = self.alpha * self.dLogEdalpha(); 
        end
        
        function res = dLogEdalpha(self)
            res = -(self.H(self.s) - self.ApproximateGalpha(self.s)/(2*self.alpha));
        end
        
        function X = CrossEntropy(self,s_)
            [ds, g] = self.GetDsdU(s_);
            X = -0.5*(g'*ds);
        end
        
        function res = ApproximateGalpha(self, s, niter)
            if nargin<3, niter = 1; end
            % currently, niter is ignored.
            % cache? If[ValueQ[storedG], storedG,
            SUM = 0;
            for k=1:niter
                %v = randn(size(s));
                if isempty(self.r)
                    self.r = randn(size(s));
                end
                v = self.r;
                SUM = SUM + v'*self.ApplyInvB(s, self.ApplyA(s, v));
            end
            res = SUM / (niter * self.alpha);
% 			storedG = SUM / (niter alpha) ] ];
        end
        
        function [ds, g] = GetDsdU(self,s_)
            if ~isempty(self.storedDsdU)
                ds = self.storedDsdU{1};
                g = self.storedDsdU{2};
            else
                g = self.delU(s_);
                ds = self.ApplyInvK(s_, -g);
                self.storedDsdU = {ds,g};
            end
        end
        
        function g = delU(self,s_)
            g = -self.beta*(self.FTx - self.ApplyFTF(s_)) + self.alpha*self.delH(s_);
        end
        
        function res = ApplyFTF(self,y_)
            res = self.ApplyFT(self.ApplyF(y_));
        end
        
        function res = ApplyF(self,y_)
            % sum_r F(r) Y(s-r)
            % NB: conv2 generalizes to 1 dim
            % Also, convolution(F,Y) = correlation(reverse(F),Y)
            %res = conv2(self.f,y_);
            f_ = reshape(self.f,self.fr,self.fc);
            y_ = reshape(y_,self.xr,self.xc);
            if 0
                res = filter2(f_(end:-1:1,end:-1:1),y_);
            else
                res = conv2(f_,y_);
                res = res(1:self.xr,1:self.xc);
            end
            res = res(:);
            %ListConvolve[f, y, {1, 1}, 0];
        end
        
        function res = ApplyFT(self,y_)
            % sum_r F(r) Y(s+r)
            % correlate using filter2
            f_ = reshape(self.f,self.fr,self.fc);
            y_ = reshape(y_,self.xr,self.xc);
            if 0
                res = filter2(f_,y_);
            else
                res = conv2(f_(end:-1:1,end:-1:1),y_);
                res = res(self.fr:end,self.fc:end);
            end
            res = res(:);
            %ListCorrelate[f, y, {1, 1}, 0];
        end
        
        function u = U(self,s)
            u = 0.5 * self.beta * norm(self.x - self.ApplyF(s))^2 + self.alpha * self.H(s);
        end

        function h = H(self,ss)
            switch self.regMethod(1)
                case 1
                    M = self.regMethod(2);
                    h = sum(ss .* log(ss ./ M) - ss + M) ;
            end
        end
        
        function dh = delH(self,ss)
            switch self.regMethod(1)
                case 1
                    M = self.regMethod(2);
                    dh = log(ss/M);
            end
        end
        
        function ds = ApplyInvK(self,s_,b_)
            ds = self.ApplySqrtM(s_,...
                    self.ApplyInvB(s_,...
                        self.ApplySqrtM(s_,b_)))/self.alpha;
        end
        
        function res = ApplyM(self,ss_, xx_)
            switch self.regMethod(1)
                case 1
                    res = xx_ .* ss_;
            end
        end
        
        function res = ApplySqrtM(self,ss_, xx_)
            switch self.regMethod(1)
                case 1
                    res = xx_ .* sqrt(ss_);
            end
        end
        
        function res = ApplyInvM(self,ss_, xx_)
            switch self.regMethod(1)
                case 1
                    res = xx_ ./ ss_;
            end
        end
        
        function res = ApplyInvB(self,s_, b_)
            Afun = @(x) self.ApplyA(s_,x);
            res = self.alpha * approximateInversion.invert(Afun, b_, self.alpha, 0.05, 2000);
        end
        
        function res = ApplyA(self,s_, x_)
            res = self.ApplySqrtM(s_, self.beta*self.ApplyFTF(self.ApplySqrtM(s_,x_)));
        end
        
    end
end