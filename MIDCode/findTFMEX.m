function [ratefunc, TFSpan, pxt, px] = findTFMEX(nfilt, in_vec, numSpikes, stim, spikes, BINS, nvec, varargin)
%fprintf('In TF\n');

global HCURRENT;
if nargin > 7
    HCURRENT=varargin{1};
elseif isscalar(HCURRENT) == 0
    HCURRENT=0;
    fprintf('\nPost spike current not explicitly declared, assuming it does not exist\n');
end

%size(vec)
%Ensure the integers are actually integers...

if ~iscell(in_vec)
    if nvec>1
        for kk=1:nvec
            vec{kk}=in_vec(:,kk)';
        end
    else
        vec=cell(1,1);
        vec{1}=in_vec;
    end
else
    vec=in_vec;
end

tdim = length(vec{1});
nfilt=round(nfilt);
BINS=round(BINS);
nvec=round(nvec);
stimlength=length(stim);
tvec=cell(1,nvec);
stp=cell(1,nvec);
xmin=cell(1,nvec);
step=cell(1,nvec);
xmax=cell(1,nvec);

for zz=1:nvec
	if(norm(vec{zz}) ~= 0)
		vec{zz}=vec{zz}/(norm(vec{zz}));
	end
end

mini = -1.5*max(abs(stim));
maxi = 1.5*max(abs(stim));
stepsize=((maxi-mini)/BINS);

%TFSpan=cell(1,nvec);
%filt=zeros(1,nfilt(end));
%HCURRENT=1;  % Slap-fix for issues passing global variables through MEX files, need to address
if (HCURRENT == 0 || nvec==1)
    %fprintf('#');
    nspk=sum(spikes);
    spkt=find(spikes);
    spkt=spkt-tdim;
    %xmin = -10; %min(stp_o); 
    %xmax = 10;  %max(stp_o);
    for k=1:nvec
        % Find the stim projections
        % reverse the filter in each dimension
        tvec{k} = vec{k}(end:-1:1);
        % Compute the projection
        stp{k} = conv(stim,tvec{k},'valid');  %real(ifft(fft(stim, nn).*fft(vec{k},nn))); %[zeros(1,nfilt-1), conv(stim,tvec{k},'valid')]; %sameconv(stim', vec{k}')';
        ntrials=length(stp{k});
    end

    for kk=1:nvec
        xmin{kk}=mini;  %min(stp{kk});
        xmax{kk}=maxi;  %max(stp{kk});
        step{kk}= (xmax{kk}-xmin{kk})/BINS; %stepsize;
    end
    
    Ntrials=length(stp{1});
    %step=(xmax-xmin)/BINS;
    TFSpan=linspace(mini,maxi,BINS);

    %Build the histogram via MEX file
    [px, pxt]=getDist(nspk, Ntrials, xmin, step, BINS, nfilt, stp, spikes, nvec);
    ratefunc = zeros(1,length(px));
    for k=1:length(px)
        if px(k) > 0
            ratefunc(k) = pxt(k)/(px(k));
        end
    end
    %save ./debug-out.mat;
    if .9 > sum(pxt) || sum(pxt) > 1.1
        fprintf('\n\nsum of pxt: %f\n',sum(pxt));
        error('findTF:argCheck','\n P(x|spike) does not sum to one\n');
    end
    %fprintf('out TF\n');
    %fprintf(1,'\n exit tf\n ');
    
else
    % Add the post-spike current to the stim after each spike
   
    for k=1:(nvec-1)
        % Find the stim projections
        %vec{k}=vec{k}/norm(vec{k}); %Normalize the filters        
        % reverse the filter in each dimension
        tvec{k} = vec{k}(end:-1:1);
		stp{k} = conv(stim,tvec{k},'valid')'; 
		
		% Set up the histogram measurements
        xmin{k}= mini;  %min(stp{k});
        xmax{k}= maxi;  %max(stp{k});
        step{k}= ((xmax{k}-xmin{k})/BINS); %stepsize;       
    end

    % Generate convolution and binning parameters for post-spike current.
    % Parameters are adaptive to account for lack of normalization in 
    % the post-spike current.

    stp{nvec}=conv(spikes,vec{nvec},'valid');
    xmin{nvec}=mini;  %min(stp{nvec});
    xmax{nvec}= maxi;  %max(stp{nvec});
    step{nvec}= ((xmax{nvec}-xmin{nvec})/BINS);
        
    % Errors 
    if isnan(xmin{k})>0 || isnan(xmax{k})>0
        error('findTF:argCheck', '\n Error, xMin and/or xMax is NaN');
    end
    
    if BINS <= 0
        error('findTF:argCheck','\n Error, number of bins is zero\n');
    end

    TFSpan=linspace(mini,maxi,BINS);
    %step = TFSpan(2)-TFSpan(1); 
    Ntrials=length(stp{1});

    [px, pxt]=getDist(numSpikes, Ntrials, xmin, step, BINS, nfilt, stp, spikes, nvec);

    %end
    ratefunc = zeros(1,length(px));
    for k=1:length(px)
        if px(k) > 0
            ratefunc(k) = pxt(k)/(px(k));
        end
    end

    if .9 > sum(sum(pxt)) || sum(sum(pxt)) > 1.1
        save dump px pxt vec;
        fprintf('\n\nsum of pxt: %f\n',sum(pxt));
        error('findTF:argCheck','\n P(x|spike) does not sum to one\n');
    end

    clear stp;

end

