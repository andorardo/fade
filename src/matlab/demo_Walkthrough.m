function [F2] = demo_Walkthrough
%Demonstation of the use of Fast Automated Deconvolution of EPSCs
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

%% Load data

scriptDir = mfilename('fullpath');
scriptDir = scriptDir(1:end-length(mfilename)-1);
load([scriptDir filesep '..' filesep 'resources' filesep 'demodata.dat']);


%% Construct EPSC deconvolution object
F = FADE(20e3,demodata(1:1e4));
F.threshold = 10;
F.smoothing = 0;


%% Iterate a few times

for k=1:6
    F.update;
end


%% Plot some statistics of the larger dataset

F2 = FADE(20e3,demodata,F.f);
F2.hFig = F.hFig;
F2.threshold = 10;
F2.smoothing = 0;
F2.update;
% Add polish:
F2.update;

figure;
subplot(2,2,1);
hist(F2.amplitudes,25);
subplot(2,2,2);
dt = diff(F2.times);
hist(-exp(-dt/mean(dt)),25);
subplot(2,2,3);
normplot(F2.residual);


%% Compare with previously calculated 

% load([scriptDir filesep '..' filesep 'resources' filesep 'demo_bestf.dat']);


end