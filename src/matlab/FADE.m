classdef FADE < MEMDeconvolution
%Fast automated deconvolution of EPSCs.
%
%   Implements blind deconvolution of EPSC-type problems by alternately
%   updating signal and point spread function ('filter').
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
        samplingRate;
        data;
        baseline;
        hFig;
        doPlot = true;
        
        alphaPaths = {};
        xMeans = [];
        baselineStats = [];
        filterEstimates = [];
        residualMeans = [];
    end
    
    methods (Access=public)
        
        function obj = FADE(srate_,data_,f_,varargin)
            if nargin<3 || isempty(f_)
                f_ = FADE.filterFunc(srate_,(0:174)');
            elseif iscell(f_)
                f_ = FADE.filterFunc(srate_,(0:f_{1}-1)',f_{2:end});
            end
            x = data_ - median(data_);
            % Flip data if negative currents
            if sum(x.*sign(x)) > sum(x.*sign(-x))
                x = -x;
            end
            data_ = x;
            bline = -gaussianDetrend(x,srate_,1,true);
            % Better fit after some fluctuations removed:
            if ~isempty(which('mle')) % check to see if Stats Toolbox function 'mle' is available
                % XXXX subsample x for very large data
                % XX possibly use own pdf function
                [phat] = mle(x + bline,'Distribution','tlocationscale');
                bline = bline - phat(1);
            end
            obj = obj@MEMDeconvolution(x + bline, f_, varargin{:});

            obj.samplingRate = srate_;
            obj.data = data_;
            obj.baseline = bline;

            obj.appendStats;
        end
        
        function updatePlot(self)
            if ~self.doPlot, return; end
            if isempty(self.hFig), self.hFig = figure; end
            
            tx = (0:length(self.x)-1)/(self.samplingRate/1e3);
            tf = (0:size(self.f,1)-1)/(self.samplingRate/1e3);
            ncol = max(6,size(self.filterEstimates,2));

            figure(self.hFig);
            
            events = self.times;
            if ~isempty(events)
                dfixes = {1:2,3};
                d = diff(events);
                w = min(length(d),4);
                sz = filter(ones(w,1),1,d);
                sz = sz(w:end);
                [~,IX] = sort(sz);
                IX = IX(ceil(length(IX)/2));
                e1 = events(IX);
                e2 = e1 + max(length(self.f),events(IX+w-1)-e1);
                de = floor((e2-e1)/8);
                ixes = {1:length(tx),e1-de:e2-de};
            else
                dfixes = {1:3};
                ixes = {1:length(tx)};
            end
            for m=1:length(dfixes)
                subplot(3,3,dfixes{m}); cla;
                ix = ixes{m};
                stairs(tx(ix),self.x(ix),'k'); hold on;
                stairs(tx(ix),self.sparseS(ix),'r');
                stairs(tx(ix),self.reconstruction(ix),'m');
                title('Data and reconstructions');
                xlabel('Time (msec)');
                axis tight;
            end

            subplot(3,3,4:5); cla;
            set(gca,'NextPlot','replacechildren')
            set(gca,'ColorOrder',flipud(copper(ncol)));
            stairs(tf,self.filterEstimates);
            title('Filter estimates');
            xlabel('Time (msec)');
            axis tight;

            if size(self.filterEstimates,2)>1
                dfs = diff(self.filterEstimates,1,2);
                ssq = sum(dfs.^2,1);
                subplot(3,3,6); cla;
                semilogy(sqrt([nan ssq]),'-*');
                title('Norm of filter changes');
                xlabel('Iter');
            end

            
            subplot(3,3,7); cla;
            plot(self.residualMeans,'-*');
            hold on;
            plot(zeros(size(self.residualMeans)),'--');
            title('Mean of residuals');
            xlabel('Iter');
            
            subplot(3,3,8); cla;
            plot(self.baselineStats(1,:),'-*');
            title('Mean of baseline adjustment');
            xlabel('Iter');

            subplot(3,3,9); cla;
            plot(self.baselineStats(2,:),'-*');
            title('Std of baseline adjustment');
            xlabel('Iter');

            drawnow;
        end
        
        function update(self)
            % Perform one cycle of updates:
            % 1. Deconvulution -- estimate x given f
            % 2. Filter estimation -- estimate f given x

            % Deconvolution
            self.resetS;
            self.alpha = 1e3;
            tic; self.MaximizeEvidenceAlphaOnly; toc;

            self.updatePlot;
            
            % Baseline adjustment
            smoothResid = gaussianDetrend(self.residual,self.samplingRate,self.threshold,true);
            self.baseline = self.baseline - smoothResid;
            
            % Filter estimation
            fnew = self.filterEstimate('method','lp','refinepeak',false,'peakposition',25);
            fnew = fnew.*double(fnew>0);
            self.f = fnew;

            % Some record keeping
            self.appendStats;
            
            % Update x
            self.x = self.data + self.baseline;
            
            % Ready for the next update.
        end
        
    end

    methods (Access=private)
        
        function appendStats(obj)
            obj.xMeans(end+1) = mean(obj.x);
            obj.baselineStats(1:2,end+1) = [mean(obj.baseline);std(obj.baseline)];
            obj.filterEstimates(:,end+1) = obj.f;
            obj.residualMeans(end+1) = mean(obj.residual);
        end
        
    end
    
    methods (Static)
        
        function f = filterFunc(srate,t,tau,tau1)
            % srate is the sampling rate
            scale = srate/20e3;
            if nargin<3, tau=scale*9.8; end
            if nargin<4, tau1=scale*6.9; end
            et1 = exp(-t./tau1);
            et = exp(-t./tau);
            tt = tau + tau1;
            f = ((1 - et1) .* et) * (tt/tau1)^(tau1/tau)*tt/tau;
        end
        
    end

    
    
    
end