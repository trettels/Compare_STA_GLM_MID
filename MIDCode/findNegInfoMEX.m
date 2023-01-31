%% Compute the Negative Information: Used in the simulated annealing
% algorithm 
function[neginfo]=findNegInfoMEX(stim,spikes,inVec, nfilt, BINS,nvec, varargin)
global HCURRENT;
if nargin > 6
    HCURRENT=varargin{1};
elseif isscalar(HCURRENT) == 0
    HCURRENT=0;
    fprintf('\nPost spike current not explicitly declared, assuming it does not exist\n');
end

% Make sure the vector is in the proper format:
if iscell(inVec)
    vec=inVec;
elseif nvec==1
    vec{1}=inVec;
else
    vec=cell(1,nvec);
    [m,n]=size(inVec);
    if nvec ~=1 && m<n
        inVec=inVec';
        fprintf('input vector improperly aligned, taking the transpose...');
    end
    for k=1:nvec
        vec{k}=inVec(:,k)';
    end
end

nvec=round(nvec);
%Ensure the integers are actually integers...
nfilt=round(nfilt);
BINS=round(BINS);

% Get the Probability distributions...
[~,~, pxt, px] = findTFMEX(nfilt, vec, sum(spikes), stim, spikes, BINS,nvec);

    % Compute the negative information
[neginfo]=negInfo(px, pxt, length(pxt));
clear px pxt;

