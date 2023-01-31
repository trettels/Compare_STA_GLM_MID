function [jackfilt, partjackfilt, midrf, dom] = maxInform(stimtemp, spikes, nfiltprime, BINS, stafilt, startv, options,nvec)
% The primary function for the matlab implementation of the MID
% algorithm.  
DSF=options.DSF;
global HCURRENT;
dtime=zeros(1,options.MaxRep);
nfilt = nfiltprime;
stim=stimtemp;
resamp = options.Resamp;
%% Pre-Processing
% Downsample, if necessary
% Downsampling of the spike train, without loosing spikes
if DSF ~=1
    nfilt = floor(nfilt/DSF);
    t=find(spikes);
    st=downsample(stim, DSF);
    sp=zeros(size(st)); %temp spike-train vector shell
    for i=1:length(t)
        if floor(t(i)/DSF) < 1
            sp(1)=sp(1)+1;
        else
            sp(floor(t(i)/DSF))=sp(floor(t(i)/DSF))+1; % Bin the spikes
        end
    end
else
	st=stim;
    sp=spikes;
end

% Assign the size variables
stimlength=length(st);

% Setup the seed-vector
bestv=cell(1,nvec);

[m, n] = size(startv);
if m<n
	startv=startv';
	fprintf('\nRe-orienting starting vector\n');
end
for p=1:nvec    
    bestv{p}= startv(:,p)';
end
testv=bestv;
initv=testv;

info_b = findNegInfoMEX(st,sp,bestv,nfilt,BINS,nvec);
fprintf(1,'Starting Info = %7.5f\n',-info_b);
ndim=ndims(st);
%% Mex Jackknife
% Set up for a jackknife-like run, splitting the stimulus and
% spikes up into 8 parts and running 8 parallel annealing runs,
% each running over 7 of the 8 (recombined) fragments
fin=99;
while (fin > 0) && (resamp > 1)
	sti = cell(1,resamp);
	spik=cell(1,resamp);
	for i=0:resamp-1
		s=i*(floor(stimlength/resamp));
		sti{i+1} = st(s+1:s+(floor(stimlength/resamp)));   
		spik{i+1} = sp(s+1:s+(floor(stimlength/resamp)));
	end	
	iter_temp=resamp;
	temp_spk=spik;
	for i=1:iter_temp
		temp_spk{i}(1:nfilt)=0;
		if(sum(temp_spk{i}) <=10)
			resamp=resamp-1;
			fin=99;
			break;
		else
			fin=-99;
		end
	end
end
tic
partjackfilt = cell(1,resamp);
for jack=1:resamp
    stimTrain = [];
    spikeTrain=[];
    spikeTest=[];
    stimTest=[];
    % Set up the train and test sets...
    for i=1:resamp
        if i~=jack && resamp > 1
            stimTrain = cat(ndim,stimTrain, sti{i}); % Concatinates the stimulus segments along the time-dimension
            temp = spik{i};
            if i==jack+1
	            temp(1:nfilt)=0;
	        end
            spikeTrain =[spikeTrain, temp];
        end
        if i==jack && resamp > 1
            stimTest = sti{i};
            spikeTest = spik{i};
            spikeTest(1:nfilt)=0;
        elseif jack==0 || resamp == 1
            stimTrain=st;
            spikeTrain=sp;
            stimTest = st;
            spikeTest=sp;
        end
    end
   
    fprintf(1,'\n Starting set %u of %u, number of spikes: %u (train), %u (test) \n',jack, resamp, sum(spikeTrain), sum(spikeTest));

    % Initialize the vectors and comparison cases
    q=initv;
    bestq=q;
    bestq_pass=bestq;
    fprintf('Finding information of train set\n');
    infob_train = findNegInfoMEX(stimTrain,spikeTrain,initv,nfilt,BINS,nvec);
    infoq_train = infob_train;
    infob_pass=infob_train;
    infoq_pass=infoq_train;
    fprintf('finding information of test set\n');
    
    save debug-out;
    infob = findNegInfoMEX(stimTest,spikeTest,bestq,nfilt,BINS,nvec);
    fprintf('\n starting annealing\n');


    
    %% Low-res Search
    y=0;
    T = options.InitTemp;
    minTemp=options.FinalTemp;
    while (y<options.MaxRep) && (T > options.FinalTemp)

        ticky=tic;
        y=y+1;
        [q, T] = testanneal(options.CoolSched, T, options.MaxConsRej, 20, options.MaxTries, options.MaxIter, minTemp, options.Heat, -inf, q, stimTrain, stimTest, spikeTrain, spikeTest, nfilt, BINS, infob_pass, infoq_pass, options.InitTemp, bestq_pass, prod(nfilt),nvec);
		
        info_q = findNegInfoMEX(stimTest,spikeTest,q,nfilt,BINS, nvec); %get the negative info for the normed new vector
        infoq_train = findNegInfoMEX(stimTrain,spikeTrain,q,nfilt,BINS, nvec);
        
        %Check to see if the new vector (normed) is better than the prior best
        if  info_q < infob
            bestq=q;
            bestq_pass=bestq;
            infob=info_q;
            infob_train = infoq_train;
            infob_pass=infob_train;
            infoq_pass=info_q;
            infot=findNegInfoMEX(st,sp,bestq,nfilt,BINS,nvec);
            fprintf(1,'\n\t\t (%u.%u) >> Updated, new info = %7.5f\n', jack, y, -infot);
        else
            qtemp=1-2*rand(nfilt,nvec);
            
            for kk=1:nvec
                q{kk}=(q{kk} +(1+T)*qtemp(:,kk)');
            end
            infoq_pass=findNegInfoMEX(stimTrain,spikeTrain,q,nfilt,BINS, nvec);
            infob_pass=infoq_pass;
            bestq_pass=q;
            fprintf(1,'\n (%u.%u)\t>>> Current Temperature: %7.7f\n',jack, y, T);
        end
        tocky=toc(ticky);
        dtime(y)=tocky/60;
        fprintf('\n\nIteration Run Time: %7.7f minutes\n\n',dtime(y));
    end
    for zz=1:nvec
	    partjackfilt{jack}{zz}=interp(bestq{zz}, DSF);  % Output the best vector found
	end
	
end

%% Ensure vectors are oriented properly and that all similar vectors are in the same vector position
for zz=1:nvec-1
    sign1=sign(dot(stafilt/norm(stafilt), partjackfilt{1}{zz}/norm(partjackfilt{1}{zz})));
    % Ensure all the vectors are going the same way
    for k=2:resamp
        stest1=sign(dot(partjackfilt{1}{zz}/norm(partjackfilt{1}{zz}), partjackfilt{k}{zz}/norm(partjackfilt{k}{zz})));
        if sign1 ~= stest1
            fprintf('Flipping vector %i in sample %i\n',zz,k);
            partjackfilt{k}{zz}=-partjackfilt{k}{zz};
        end
    end
end

   
temp=cell(1,nvec);
jackfilt=cell(1,nvec);

for k=1:nvec
    temp{k}=0;
    for r=1:resamp
        part=partjackfilt{r};
        temp{k} = temp{k}+part{k};
    end
    jackfilt{k}=temp{k}/resamp;
end

% Orient the same direction as STA
for zz=1:nvec-HCURRENT
    if sign( dot( stafilt/norm(stafilt), jackfilt{zz}/norm(jackfilt{zz}) ) ) ~= sign( dot( stafilt/norm(stafilt), stafilt/norm(stafilt) ) )
        jackfilt{zz}=-jackfilt{zz};
        fprintf('\nFlipping MID Filter\n');
    end
end
toc

[midrf, dom, ~, ~] = findTFMEX(nfiltprime, jackfilt, sum(spikes), stimtemp, spikes, BINS, nvec);
if nvec > 1
    dims=[BINS];
    for dd=2:nvec
        dims=[dims, BINS];
    end
    midrf=reshape(midrf,dims);
end