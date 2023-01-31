function [ratefunc, STA, STCFilt, eVals, eVec, STC, LAMBDA] = genSTASTC(stimulus, spike, duration, BINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenSTASTC:
%    Generate the STA and STC functions for a model
%    Sean Trettel, Fairhall Lab, Biophysics Dept, University of Washington
%
% This function is to generate and return the Spike Triggered Average (STA)
% and the Spike Triggered Covariance (STC), the corresponding
% eigenvalues and first eigenvector for a given stimulus vector (stim) and
% spike-train vector (spikes), the normal vector to each stimulus block,
% the matrix of stimulus blocks, the spike times and finally the triggering
% stimulus blocks.  Duration is the size of the stimulus vector (Length of
% the trial) and t_dim is the time-step size.
%
% The spikes vector should be a boolean vector with 'spikes(i)=1'
% indicating a spike at time 'i'.
%
% Modified on 03 JUN 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants, Paramaters and Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the window size of the spike triggering blocks (in ms)
blockSize = duration;  
% Drop the spikes with insufficient history by setting all timepoints prior
% to (blocksize) ms to zero
spike(1:blockSize) = 0;
%testSpikes = spikes;
% Extract the triggered spike times and count them
if find(spike>1) > 0
	n=1;
	for i=1:length(spike)
		j=1;
		while j <= spike(i)
			spikeTimes(n)=i;
			n=n+1;
			j=j+1;
		end
	end
else
	spikeTimes = find(spike);
end
numSpikes = sum(spike);
% Intitialize the matrix of triggering stimuli vectors
trigStim = zeros(numSpikes,blockSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the STA and recover the nonlinear filter, and reorganize
% stimulus vector into a matrix of(blocksize) millisecond sample stimulus
% row vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim = zeros(length(stimulus)-blockSize,blockSize);
% Generate the set of spike generating stimulus blocks, and split the
% stimulus vector into a matrix of stimulus blocks of size (blockSize)


for k = 1:(length(stimulus)-blockSize)
    stim(k,:) = stimulus(1,k:(k+blockSize)-1);
end

for k = 1:numSpikes
    trigStim(k,:) = stimulus(1,(spikeTimes(k)-blockSize):(spikeTimes(k)-1));
end

% Find the STA vector
STA = mean(trigStim-mean(stimulus)); %-mean(stimulus));
STA = STA; %./norm(STA);

fprintf('\n found STA, starting STC\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the STC and conduct eigenvalue and eigenvector analysis on it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the covariance of the stimuli and of the triggering stimuli
LAMBDA=cov(trigStim);
STC = LAMBDA - cov(stim);
% Pull out the Eigenvalues and Eigenvectors
[eVec, eigVal] = eig(STC);
eVals = diag(eigVal);
STCFilt=eVec(:,find(abs(eVals)==max(abs(eVals))));
STCFilt=STCFilt/norm(STCFilt);

fprintf('\n STC and corresponding eigenvalue analysis done, starting i/o functions\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover the Nonlinearity (currently a direct fitted copy from old code):
% Need to re-scale so that it centers properly.  Way to do it
% non-artificially?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratefunc=cell(1,2);
for i=1:2
    if i==1        
        vec=STA;
    else
        vec=STCFilt';
    end  
    vec = vec./norm(vec);
    cvec=vec(end:-1:1);
    stimproj=conv(stimulus,cvec,'valid')';
    stsproj=zeros(1,numSpikes);
	    for k=1:numSpikes 
        stsproj(k)= dot(trigStim(k,:), vec)/dot(vec,vec);% P(x|Spike)      
    end
    mini=-1.5*max(abs(stimulus));
    maxi=1.5*max(abs(stimulus));
    [nstim,x]=hist(stimproj,linspace(mini,maxi,BINS));    
    [nsts,x]=hist(stsproj,linspace(mini,maxi,BINS));
    ratefunc{i}=nsts./(nstim+eps); % divide histograms to get rate
end