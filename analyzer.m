function [sta, mid, glm, starttime, endtime] = analyzer(stim, tsp, varargin)
%% Analyzer function:  Takes stimulus and spike times with optional options structure and returns 
% MID, GLM and STA data structures with filters, rate functions and conditional intensity functions
starttime=datevec(now);
if ~isempty(varargin)
    options=varargin{1};
else
    options = struct('MaxTries',2);
end

%% Set up the options structure

% MaxTries is used in the simulated annealing algorithm, and sets the number of iterations run
% before decreasing the temperature  Default value is 2.
if ~isfield(options,'MaxTries') || isempty(options.MaxTries)
    options.MaxTries=2;
end

% MaxConsRej is used in the simulated annealing algorithm, and is the number of consecutive
% rejected iterations allowed before reducing the heating coefficient and reseting the vector to
% the best vector found in the current annealing run.  Default is 30.
if ~isfield(options, 'MaxConsRej') || isempty(options.MaxConsRej)
    options.MaxConsRej=30;
end

% InitTemp is used in the simulated annealing algorithm and is the inital temperature used.
% Default is 0.05.
if ~isfield(options, 'InitTemp') || isempty(options.InitTemp)
    options.InitTemp=0.05;
end

% MaxIter is used in the simulated annealing algorithm, and is the maximum number of anealing 
% iterations within a single annealing run. Defualt is 250 iterations.
if ~isfield(options,'MaxIter') || isempty(options.MaxIter)
    options.MaxIter=250;
end

% CoolSched is used in the simulated annealing algorithm, and is the scalar by which the temperature
% is multiplied when the temperature is cooled.  Default is 0.98.
if ~isfield(options,'CoolSched') || isempty(options.CoolSched)
    options.CoolSched=.98;
end

%Heat is used in the simulated annealing algorithm, and is the inital value for the heating
% coefficient.  The heating coefficient is used to increase the temperature when jittering the test
% vector when the vector is not moved during a minimization.  Default is 1.2.
if ~isfield(options, 'Heat') || isempty(options.Heat)
    options.Heat= 1.2;
end

% FinalTemp is used in the simulated annealing algorithm and is the temperature at which no more
% annealing runs will occur, regardless of MaxRep.  Default is 5e-12.
if ~isfield(options, 'FinalTemp') || isempty(options.FinalTemp)
    options.FinalTemp=5e-12;
end

% MaxRep is used in the simulated annealing algorithm, and is the maximum number of simulated
% annealing runs allowed to occur, regardless of the temperature.  Default is 10.
if ~isfield(options, 'MaxRep') || isempty(options.MaxRep)
    options.MaxRep=10;
end

% resamp is used in the simulated annealing algorithm, and is the number of times to resample the
% data for jackknife runs.  Default is 1 (no resampling)
if ~isfield(options, 'Resamp') || isempty(options.Resamp)
    options.Resamp=1;
end

% DSF is used in the simulated annealing algorithm as the scale by which to downsample the stimulus
% and spike train.  Output is resampled prior to passing back to 'analyzer', and will be on the 
% same scale as the original input. Default is 1 (No downsampling)
if ~isfield(options, 'DSF') || isempty(options.DSF)
    options.DSF=1;
end

% klength is used in the STA, simulated annealing and GLM algorithms, and is the length in frames 
% of the filters being searched for.  Default is 100 frames.
if ~isfield(options, 'KLength') || isempty(options.KLength)
    options.KLength=100;
end

% nfilts is the number of filter vectors to search for in the simulated annealing and spike
% triggered covariance algorithms.  Default is 1.
if ~isfield(options, 'NFilts') || isempty(options.NFilts)
    options.NFilts=1;
end

% MIDMultiRun is a flag which controls whether the MID algorithm is run either once at full
% dimensionality (N=NFilts) (MIDMultiRun < 1) or multiple times, once each for NFilts >= N > 0.
% Default is -1 (only run for N=NFilts)
if ~isfield(options, 'MIDMultiRun') || isemtpy(options.MIDMultiRun)
    options.MIDMultiRun = options.NFilts
end

% Bins is the number of bins per filter used in creating probability distributions in the MID
% algorithm.  Default is 25 bins.
if ~isfield(options, 'Bins') || isempty(options.Bins)
    options.Bins=25;
end

% FirstPeak is used in the GLM algorithm, and defines the location in frames of the first basis
% peak.  Default is 1.
if ~isfield(options, 'KFirstPeak') || isempty(options.KFirstPeak)
    options.KFirstPeak=1;
end

% Last peak is used in the GLM algorithm and defines the location in frames of the last basis peak.
% Default is one half of klength, rounded.
if ~isfield(options,'KLastPeak') || isempty(options.KLastPeak)
    options.KLastPeak=round(0.5*options.KLength);
end

% Linearity is used in the GLM algorithm and determines how linearly the basis peaks are distributed
% between the first and last peaks.  A higher number gives a more linear dsitribution.  Default is 5.
if ~isfield(options,'KLinearity') || isempty(options.KLinearity)
    options.KLinearity=5;
end

% KBasis is used in the GLM algorithm, and is the number of raised cosine basis functions to use in
% creating the stimulus kernel.  Default is 10.
if ~isfield(options,'KBasis') || isempty(options.KBasis)
    options.KBasis=10;
end

% FirstPeak is used in the GLM algorithm, and defines the location in frames of the first basis
% peak of the post-spike current.  Default is 1.
if ~isfield(options, 'HFirstPeak') || isempty(options.HFirstPeak)
    options.HFirstPeak=1;
end

% Last peak is used in the GLM algorithm and defines the location in frames of the last basis peak
% of the post-spike current.  Default is one half of klength, rounded.
if ~isfield(options,'HLastPeak') || isempty(options.HLastPeak)
    options.HLastPeak=round(0.5*options.KLength);
end

% HLinearity is used in the GLM algorithm and determines how linearly the post-spike basis peaks are
% distributed between the first and last peaks.  A higher number gives a more linear dsitribution.
% Default is 5.
if ~isfield(options,'HLinearity') || isempty(options.HLinearity)
    options.HLinearity=5;
end

% HNBasis is used in the GLM algorithm, and is the number of raised cosine basis functions to use in
% creating the post-spike current kernel.  Default is 10.
if ~isfield(options,'HBasis') || isempty(options.HBasis)
    options.HBasis=10;
end

% AbsRef is the duration of absolute refractory period assumed in frames.  Default is 1 frame.
if ~isfield(options,'AbsRef') || isempty(options.AbsRef)
    options.AbsRef=1;
end

% RefreshRate is used in the GLM algorithm to modulate the spike rate during the maximum likelihood
% estimation process.  (RefreshRate / DTSim) should be of the same order of magnitude as the stimulus
% sampling rate.  Default is 100.
if ~isfield(options, 'RefreshRate') || isempty(options.RefreshRate)
    options.RefreshRate=100;
end

% DTSim is used in the GLM algorithm in order to subdivide stimulus frames when calculating spike
% times in the maximum likelyhood estimation. (RefreshRate / DTSim) should be of the same order of 
% magnitude as the stimulus sampling rate.  Default is 0.1.
if ~isfield(options, 'DTSim') || isempty(options.DTSim)
    options.DTSim=0.1;
end


% GLMIter is the maximum number of iterations allowed in the maximum liklihood optimization
% algorithm.  Default value is 500.
if ~isfield(options, 'GLMIter') || isempty(options.GLMIter)
    options.GLMIter = 500;
end

% GLMTol is the tolerance at which the maximum likelihood function will consider itself complete.
% Default value is 5e-5.
if ~isfield(options, 'GLMTol') || isempty(options.GLMTol)
    options.GLMTol = 5e-5;
end

% Hcurrent is used in the MID algorithm, it acts as the toggle as to whether spike history should be
% accounted for.  HCurrent < 1 indicates that no spike history component should be considered, 
% and HCurrent >= 1 indicates that a spike history component should be searched for as the last
% vector in a cell-array of filters. Default is 0;
if ~isfield(options,'HCurrent') || isempty(options.HCurrent)
    options.HCurrent=0;
end

%%Prep mex files, paths and global variables

%Initialize paths
setpaths;

%Compile mex code if necessary
initialize_mexcode;
if ~exist('getDist')
    mex ./MIDCode/mexcode/getDist.c -outdir ./MIDCode/mexcode;
end
if ~exist('linearMinimizer')
    mex ./MIDCode/mexcode/linearMinimizer.c -outdir ./MIDCode/mexcode;
end
if ~exist('negInfo')
    mex ./MIDCode/mexcode/negInfo.c -outdir ./MIDCode/mexcode;
end
if ~exist('testanneal')
    mex ./MIDCode/mexcode/testanneal.c -outdir ./MIDCode/mexcode;
end

%Set up global variables
global HCURRENT;
HCURRENT=options.HCurrent;
global RefreshRate;
RefreshRate=options.RefreshRate;
global KLASTPEAK;
global KFIRSTPEAK;
global KLINEARITY;
global KBASIS;
KLASTPEAK=options.KLastPeak;
KFIRSTPEAK=options.KFirstPeak;
KLINEARITY=options.KLinearity;
KBASIS=options.KBasis;
global HLASTPEAK;
global HFIRSTPEAK;
global HLINEARITY;
global HBASIS;
global ABSREF;
HLASTPEAK=options.HLastPeak;
HFIRSTPEAK=options.HFirstPeak;
HLINEARITY=options.HLinearity;
HBASIS=options.HBasis;
ABSREF=options.AbsRef;

%% Derived stats

% Orient the stimulus properly.
[m n]=size(stim);
if(m > n)
    stim = stim';
end
if (m > 1) && (n > 1)
    error('Analyzer:argChk','\n\t\tThis code is not set up to handle multiple dimensions of stimuli. Sorry.\n');
    return;
end
slen=length(stim);
%Generate the spike train
spikes=zeros(1,length(stim));
for cc=1:length(tsp)
    index=round(tsp(cc));
    if index==0
        index=1
    elseif index > slen
        index=length(spikes);
    end
    spikes(index)=spikes(index)+1;
end

%% Generate the STA and STC data.  The resultant structures are as follows:
%
% sta.RateFunc: cell array of the nonlinear rate functions for the STA and STC (1D) filters.  Cell 1
%    is the STA, cell 2 is the first STC filter
%
% sta.Filt: The spike triggered average filter
%
% sta.STCFilt: The first spike-triggered covariance filter, uses the STCMatrix to calculate.
%
% sta.eVals: Eigenvalues of sta.STCMatrix
%
% sta.eVec: Eigenvectors of sta.STCMatrix
%
% sta.STCMatrix: The matrix determined by subtracting the covariance of the spike-triggering
%    stimulus matrix by the covariance of the full stimulus [ = cov(Lambda) - cov(stim)]
%
% sta.Lambda: The covariance of the spike-triggering stimulus without subtracting the covariance of
%    the full stimulus.
fprintf('\nBeginning STA...');
[sta.RateFunc, sta.Filt, sta.STCFilt, sta.eVals, sta.eVec, sta.STCMatrix, sta.Lambda] = genSTASTC(stim, spikes, options.KLength, options.Bins);

%% Generate the GLM fit
%
% glm.NegLogLiVal: the final negative log likelihood value found in the GLM algorithm.
%
% glm.Filt: The stimulus filter (glm.Filt{1}) and post-spike current(glm.Filt{2}) found by the GLM
%
% glm.DC: The DC offset found by the GLM algorithm
%
% glm.OrigStruct: The original structure produced by Jonothan Pillow's GLM code.
if options.GLMIter > 0
    gg0 = makeFittingStruct_GLM(sta.Filt',options.DTSim);  % projects sta into basis for fitting k
    gg0.tsp = tsp;  % Insert spikes into fitting struct
    gg0.tspi = 1;   % First spike to use (you can ask it to ignore the first "n" spikes)
    opts = {'display', 'iter', 'TolFun',options.GLMTol,'maxiter', options.GLMIter};
    fprintf('\nBeginning GLM fitting...');
    [gg, glm.NegLogLiVal] = MLfit_GLM(gg0,stim',opts); % do ML (requires optimization toolbox)
    % REpack the gg structure to mirror the sta and mid structures
    glm.Filt{1}=gg.k';

    temp=min(gg.iht(end), options.KLength);
    lng=length(gg.iht(1:find(gg.iht>=temp,1)));
    glm.Filt{2}=resample((gg.ihbas*gg.ih)',round(temp), lng, 0)';
    if length(glm.Filt{2} < options.KLength)
        lnt=options.KLength - length(glm.Filt{2});
        tail=zeros(lnt,1);
        glm.Filt{2}=[ glm.Filt{2}; tail]';
    end
    
    glm.DC=gg.dc;
    glm.OrigStruct=gg;
    [glm.MeanPerc, glm.MeanSuccess, glm.FitSpikes, glm.CondInt] = compareModel(stim, spikes, glm, options,2, 1);
    [glm.dkl]=findNegInfoMEX(stim, spikes, glm.Filt, options.KLength, options.Bins, 2, 1);
    [rftemp,~, glm.pxt, glm.px] = findTFMEX(options.KLength, glm.Filt, sum(spikes), stim, spikes, options.Bins, 2);
    glm.pxt=reshape(glm.pxt, options.Bins, options.Bins);
    glm.px=reshape(glm.px, options.Bins, options.Bins);
    glm.RateFunc=reshape(rftemp, options.Bins, options.Bins);
else
    fprintf('\nIgnoring GLM...\n');
end
%% Generate the MID Fit
%
% mid.filt{NFilts}: The MID-fitted filters. If HCurrent > 0, mid.filt{1:NFilts-1} are the stimulus
%    filters and mid.filt{NFilts} is the post-spike current
%
% mid.part{resamp}: The cell array of filters found in the resampled runs each mid.part{n} is a 
%    structure identical to mid.filt
%
% mid.RateFunc: Matrix or vector containing the rate function found during the MID algorithm.  If 
%    NFilts is greater than 2, manual reshaping is required.
%
% mid.RFDomain: The cell array of domain vectors for each dimension of the rate function.
if(length(options.MIDMultiRun) == 1)
    init1= 1-2*rand(floor(options.KLength/options.DSF),options.NFilts);
    [mid.Filt, mid.Part, mid.RateFunc, mid.RFDomain] = maxInform(stim, spikes, options.KLength, options.Bins, sta.Filt, init1, options, options.MIDMultiRun);
    mid.infograph=load('anneal.txt', '-ascii');
else
    mid=cell(1,length(options.MIDMultiRun));
    for NN=options.MIDMultiRun
        init1= 1-2*rand(floor(options.KLength/options.DSF),NN);
        [mid{NN}.filt, mid{NN}.part, mid{NN}.rf, mid{NN}.dom] = maxInform(stim, spikes, options.KLength, options.Bins, sta.Filt, init1, options, NN);
        mid{NN}.infograph=load('anneal.txt', '-ascii');
    end
end

%% Get D_kl for all three filters
[sta.dkl]=findNegInfoMEX(stim, spikes, sta.Filt, options.KLength, options.Bins, 1, -1);
[sta.STCdkl]=findNegInfoMEX(stim, spikes,sta.STCFilt, options.KLength, options.Bins, 1, -1);
if(length(options.MIDMultiRun) == 1)
    [mid.dkl]=findNegInfoMEX(stim, spikes, mid.Filt, options.KLength, options.Bins, options.MIDMultiRun, options.HCurrent);
else
    for NN=options.MIDMultiRun
        [mid{NN}.dkl]=findNegInfoMEX(stim, spikes, mid{NN}.Filt, options.KLength, options.Bins, NN, options.HCurrent);
    end
end


%% Get Distributions
%
% *.pxt: the spike conditioned distribution [ P(x|spike) ]
%
% *.px: the filtered stimulus distributions
[~,~, sta.pxt, sta.px] = findTFMEX(options.KLength, sta.Filt, sum(spikes), stim, spikes, options.Bins, 1);
[~,~, sta.STCpxt, sta.STCpx] = findTFMEX(options.KLength, sta.STCFilt, sum(spikes), stim, spikes, options.Bins, 1);

if(length(options.MIDMultiRun) == 1)
    [~,~, mid.pxt, mid.px] = findTFMEX(options.KLength, mid.Filt, sum(spikes), stim, spikes, options.Bins, options.MIDMultiRun);
    if(options.MIDMultiRun==2)
        mid.pxt=reshape(mid.pxt, options.Bins, options.Bins);
        mid.px=reshape(mid.px, options.Bins, options.Bins);
    elseif(options.MIDMultiRun > 2)
        fprintf('\nMore than two probability dimensions, please check P(x) and P(x|spike)\n');
        dims=[options.Bins];
        for dd=2:options.MIDMultiRun
            dims=[dims, options.Bins];
        end
        mid.pxt=reshape(mid.pxt, dims);
        mid.px = reshape(mid.px, dims);
    end
else
    for NN=options.MIDMultiRun
        [~,~, mid{NN}.pxt, mid{NN}.px] = findTFMEX(options.KLength, mid{NN}.Filt, sum(spikes), stim, spikes, options.Bins, NN);
        if(NN==2)
            mid{NN}.pxt=reshape(mid{NN}.pxt, options.Bins, options.Bins);
            mid{NN}.px=reshape(mid{NN}.px, options.Bins, options.Bins);
        elseif(options.NN > 2)
            fprintf('\nMore than two probability dimensions, please check P(x) and P(x|spike)\n');
            dims=[options.Bins];
            for dd=2:NN
                dims=[dims, options.Bins];
            end
            mid{NN}.pxt=reshape(mid{NN}.pxt, dims);
            mid{NN}.px = reshape(mid{NN}.px, dims);
        end

    end
end


%% Test the Models
%
% *.MeanPercent: The mean percentage of correct spike bins 
%
% *.MeanSuccess:  The mean number of correct spike bins
%
% *.FitSpikes: An array of the spike trains produced during the model testing
%
% *.CondInt: The conditional intensity function for the fitted model
%{
options.RFIndex=1;
[sta.MeanPerc, sta.MeanSuccess, sta.FitSpikes, sta.CondInt] = compareModel(stim, spikes, sta, options,1, 1);
options.RFIndex=2;
stcTemp=struct('Filt', sta.STCFilt, 'RateFunc',sta.RateFunc{2});
[sta.STCMeanPerc, sta.STCMeanSuccess, sta.STCFitSpikes, sta.STCCondInt] = compareModel(stim, spikes, stcTemp, options,1, 1);
clear stcTemp;

if(length(options.MIDMultiRun) == 1)
    [mid.MeanPerc, mid.MeanSuccess, mid.FitSpikes, mid.CondInt] = compareModel(stim, spikes, mid, options,options.MIDMultiRun, 1);
else
    for NN=options.MIDMultiRun
        [mid{NN}.MeanPerc, mid{NN}.MeanSuccess, mid{NN}.FitSpikes, mid{NN}.CondInt] = compareModel(stim, spikes, mid{NN}, options,NN, 1);
    end
end
%}
glm=struct([]);
endtime=datevec(now);
fprintf('\n\n\tStart time: %s\n\tEnd time: %s\n',datestr(starttime),datestr(endtime));