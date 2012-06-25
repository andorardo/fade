classdef approximateInversion
%Implementation of approximate matrix inversion developed by John Skilling
%and Mark Gibbs.
%   
%   Part of the package that accompanies the following article:
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
    
    methods (Static)
        
        function ts = test(dim,noPlot)
            if nargin<1
                dim = 5;
            end
            if nargin<2
                noPlot = false;
            end
            nelem = 2*dim;
            A = sparse(randi(dim,nelem,1),randi(dim,nelem,1),rand(nelem,1)*nelem-nelem/2,dim,dim);
            A = A+A'; A = A/2;
            % posdef!
            tEigs=tic;
            ev = eigs(A,1,'sa');
            toc(tEigs);
            A = A - speye(dim)*(ev(1)-sqrt(eps));
            u = dim*rand(dim,1);
            Theta = 0.01;
            
            ts = nan(3,1);
            
            normal_invert = nan(dim,1);
            if dim<8001
                iStart=tic;
                normal_invert = full(A + Theta*speye(dim))\u;
                ts(1)=toc(iStart);
            end
            
            full_invert = nan(dim,1);
            if dim<46001
                iStart=tic;
                full_invert = (A + Theta*speye(dim))\u;
                ts(2)=toc(iStart);
            end
            
            iStart=tic;
            approx_invert = approximateInversion.invert(@(h) A*h, u, Theta, 0.05, 1000);
            ts(3)=toc(iStart);
            
            if ~noPlot
                fprintf('Regular (full) invert took   %fs\n',ts(1));
                fprintf('Regular (sparse) invert took %fs\n',ts(2));
                fprintf('Skilling/Gibbs invert took   %fs\n',ts(3));

                figure(1); clf;
                subplot(3,1,1);
                plot(1:dim,full_invert,'b-',1:dim,approx_invert,'g--');
                legend('full','skilling');
                subplot(3,1,2);
                plot(full_invert-approx_invert);
                subplot(3,1,3);
                semilogy(abs((full_invert-approx_invert)./full_invert));

                norm(full_invert-approx_invert),
            end
        end
        
        function [dims,ts] = benchmark
            dims = unique(round(logspace(1,5.2,100)));
            m = length(dims);
            ts = nan(3,m);
            for i=1:m
                disp(['i:' num2str(i) ' dims(i):' num2str(dims(i))]);
                ts(:,i) = approximateInversion.test(dims(i),true);
            end
            
            figure(2);
            loglog(dims,ts);
        end
        
    end
    
    methods (Static)
        
        function res = invert(Afun_, u, Theta, tol, Kmax, varargin)
            [r,c] = size(u);
            Afun = @(h) Afun_( reshape(h,r,c) );
            u = u(:);
            
            g = u;
        	u2 = u'*u;
        	TMain = [];
            TOff = [];
        	TbarMain = [];
            TbarOff = [];
        	Q = 1.0;
            deltaQ = 2*tol;

            verb = 0; % should be from opts/varargin
            opt_SquareRoot = false;
            
            t0 = tic;

        	K = 1;
            Alpha(K) = u2;
        	es = g/norm(g);
            while K < Kmax && abs(deltaQ/Q) > tol
                if K == 1
                    h = g;
                else
                    h = g + Gamma(K-1) * h;
                end
                gh(K) = g'*h;
                Ah = Afun(h);
                Ah = Ah(:);
                Lambda(K) = Alpha(K) / (h'*Ah);
                hAh(K) = Alpha(K)/Lambda(K);
                g = g - Ah .* gh(K)./hAh(K);
                g2 = g'*g;
                Alpha(K+1) = g2;
                es(:,end+1) = g./sqrt(g2);
                Gamma(K) = Alpha(K+1) / Alpha(K);

                if K == 1
                    TMain(end+1) = hAh(1)/gh(1) + Theta;
                else
                    TMain(end+1) = hAh(K)/gh(K) + Gamma(K-1)*hAh(K-1)/gh(K-1) + Theta;
                    TOff(end+1) = -sqrt(Alpha(K)/Alpha(K-1)) * hAh(K-1)/gh(K-1);
                    TbarOff(end+1) = -hAh(K-1) * hAh(K)/Alpha(K-1);
                end
                TbarMain(end+1) = ((1 + Gamma(K)) * hAh(K) / Alpha(K) + Theta) * hAh(K);

                if K == 1
                    uhat = norm(u);
                    vhat = Alpha(1)/Lambda(1);
                else
                    uhat(end+1) = 0.0;
                    vhat(end+1) = 0.0;
                end
                Q       = 0.5 * uhat * approximateInversion.TriSolve(TOff,TMain,TOff, uhat);
                Qstar   = 0.5 * vhat * approximateInversion.TriSolve(TbarOff,TbarMain,TbarOff, vhat);
                deltaQ  = (0.5 * u2 - Qstar)/Theta - Q;
                if verb > 5
                    fprintf(['u.u=' num2str(u2) '; Q=' num2str(Q) '; Qstar=' num2str(Qstar) '; dQ=' num2str(deltaQ) '; dQ/Q=' num2str(deltaQ/Q) '\n']);
                end
                K = K + 1;
            end

            if opt_SquareRoot
                [TMain,TOff] = approximateInversion.CholeskyDecomposeTridiagonal(TMain,TOff);
                tmp = approximateInversion.TriSolve(TOff, TMain, 0.0 * TOff, uhat);
                if verb > 0
                    fprintf('Returning inverse sqrt.\n');
                end
            else
                tmp = approximateInversion.TriSolve(TOff, TMain, TOff, uhat);
            end
%             verb > 3 && Print["[uhat] = ", Dimensions[uhat], "; [es{1:K-1}] = ", Dimensions[Take[es,{1,K-1}]], ", [tmp] = ", Dimensions[tmp]];
            res = es(:,1:K-1) * tmp;
            res = reshape(res,r,c);
%             verb > 3 && Print["[res] = ", Dimensions[res]];
            if verb > 0 && K == Kmax
                fprintf(['Failed to converge after Kmax=' num2str(K) ' iterations; dQ/Q=' num2str(deltaQ/Q) '; ' ...
                    num2str(K/toc(t0)) ' iter/s\n']);
            elseif verb > 0
                fprintf(['Convergence after ' num2str(K-1) ' iterations; dQ/Q=' num2str(deltaQ/Q) '; ' ...
                    num2str(K/toc(t0)) ' iter/s\n']);
            end
        end
        
        function x = TriSolve(a, b, c, u)
            assert( length(b) == length(u) );
            assert( length(a) == length(c) );
            assert( length(a) == length(b)-1 );
            ulen = length(u);
            if ulen > 1
                B = zeros(ulen,3);
                B(1:end-1,1)    = a(:);
                B(:,2)          = b(:);
                B(2:end,3)      = c(:);
                d = -1:1;
                s = spdiags(B,d,ulen,ulen);
                x = s\u(:);
            else
                x = u/b(1);
            end
        end
        
        function [c, d] = CholeskyDecomposeTridiagonal(a, b)
            n = length(a);
            c = a(:); d = b(:);
            c(1) = sqrt(a(1));
            d(1) = b(1)/c(1);
            for i=2:n-1
                c(i) = sqrt(a(i)-d(i-1)^2);
                d(i) = b(i)/c(i);
            end
            c(n) = sqrt(a(n) - d(n-1)^2);
        end

    end
    
end