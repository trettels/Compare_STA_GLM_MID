function [meanperc, meansuccess, fitspikes, condint] = compareModel(stim, spikes, model, options, nvec, nrep)
%Test the efficacy of the data fitting methods

%options structure
global RefreshRate;
tdim=options.KLength;
BINS=options.Bins;
HCURRENT=options.HCurrent;
hasIH = HCURRENT;
%Model Structure
if ~iscell(model.Filt)
    filt={model.Filt};
else
    filt=model.Filt;
end
%Other variables
% Mouse = Generating Model
% Model = Fitted Model
mini=-1.5*max(abs(stim));
maxi=1.5*max(abs(stim));
step=(maxi-mini)/BINS;
scaler = 1; %;
fitspikes=cell(1,nrep);
meanperc=0;
meansuccess=0;

if hasIH > 0 && nvec > 1
	hasIH = 1;
else
	hasIH=0;
end

if ~isfield(model, 'RateFunc')
    fprintf('\nNo Nonlinear rate function specified, assuming an exponential RF\n');
    nonlin=@(xx)(exp(xx));
    dd=linspace(mini, maxi, BINS);
    rf=nonlin(dd);
    [m,n]=size(rf);
elseif iscell(model.RateFunc) && ~isfield(options, 'RFIndex')
    fprintf('\nThere are %i rate functions, which one would you like to use?\n',length(rf));
    s=input('Rate function number:');
    while ~isinteger(s)
        s=input('\nPlease enter an integer value\nRate function number:');
    end
    rf=model.RateFunc{s};
    [m,n]=size(rf);
    nonlin = @(xx) (scaler * rf( ind2sub(size(rf), xx )));
elseif iscell(model.RateFunc) && isfield(options, 'RFIndex')
    rf=model.RateFunc{options.RFIndex};
    [m,n]=size(rf);
    nonlin = @(xx) (scaler * rf( ind2sub(size(rf), xx )));
else
    rf=model.RateFunc;
    [m,n]=size(rf);
    nonlin = @(xx) (scaler * rf( ind2sub(size(rf), xx )));
end

%{
for aa=1:nrep
	dt=1;
	spkrate=1;
	stimlength=length(stim);
	model.spikes=zeros(1,stimlength);
	mouse.spikes=spikes;
	for kk=1:nvec-hasIH
		filt{kk}=filt{kk}/norm(filt{kk});
		stimp{kk} = conv(stim,filt{kk}(end:-1:1),'valid');
	end
	ntrials=length(stimp{1});
	
	%% Make the Post-spike current
	global HBASIS;
    if isempty(HBASIS)
        HBASIS=10;
        fprintf('Number of I_h basis functions to use for kernel not declared, using default of 10\n')
    end
    global HFIRSTPEAK;
    global HLASTPEAK;
    if isempty(HFIRSTPEAK)
        HFIRSTPEAK=1;
        fprintf('Position of first I_h basis function peak undeclared, using default of frame 1\n');
    end
    if isempty(HLASTPEAK)
        HLASTPEAK=round(nkt*0.33);
        fprintf('Position of last I_h basis function peak undeclared, using default of (0.33*window length)\n');
    end
    global HLINEARITY;
    if isempty(HLINEARITY)
        HLINEARITY=.5;
        fprintf('HLinearity of peak distribution undeclared, using default of .5\n');
    end
    global ABSREF;
    if isempty(ABSREF)
        ABSREF=1;
        fprintf('Absolute refractory period undeclared, using default of 1 frams\n');
    end
	ihbasprs.ncols = HBASIS;  % Number of basis vectors for post-spike kernel
    ihbasprs.hpeaks = [HFIRSTPEAK HLASTPEAK];  % Peak location for first and last vectors  [DTsim*10 2]
    ihbasprs.b = HLINEARITY;  % How nonlinear to make spacings
    ihbasprs.absref = ABSREF; %DTsim*10; % absolute refractory period 
	[iht_init,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
    if iht_init(end) > tdim
    	iht=inht_init(1:find(iht > tdim, 1));
	else
	    iht=iht_init;
    end
	if nvec > 1 && hasIH > 0
		ih=filt{nvec}'/norm(filt{nvec});
	else
		ih=zeros(size(filt{1}));
	end
	
	%% Pillow Spike Creation method (integrator)
	% ----- Set up -----
	slen=stimlength;
	VStim=cell(1,ndims(rf)-hasIH);
    for bb=1:nvec-hasIH
        Vstim{bb}=[zeros(tdim,1); stimp{bb}'; 0];
    end
    
	hlen = 1; ihhi = 0; ihthi = dt;
	nbinsPerEval = 10;  % Default number of bins to update for each spike
	minBins=5;
	nsp = 0;
	tsp = zeros(round(slen/5),1);  % allocate space for spike times
	for bb=1:nvec-hasIH
	    Vmem{bb} = interp1([0:slen+1]', Vstim{bb}, [0:dt:slen]', 'linear');
	end
	
	rlen = round(slen/dt);  % length of binned response
	
	% -------------- Compute interpolated h current ----------------------
    if HCURRENT > 0
        ih=resample(ih, length(iht), length(ih));
	    ihthi = [0:dt:iht(end)-1]';  % time points for sampling
	    ihhi = interp1(iht, ih, ihthi, 'linear', 0);
	    hlen = length(ihhi);
	end
	
	% ----- find spikes -----
	jbin = 1; % current time bin
	nsp = 0; % number of spikes
	
	tspnext = exprnd(spkrate);  % time of next spike (in rescaled time)
	rprev = 0;  % Integrated rescaled time up to current point
	while jbin <= rlen
		iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
			for zz=1:length(iinxt)
                temp = [];
                subspcdim=[];
                for bb=1:nvec-hasIH
                    temp=[temp, Vmem{bb}(iinxt(zz))];
                    subspcdim=[subspcdim, BINS];
                end
                subs=max(1, min(BINS, round((temp-mini)/step)));
                subspc=ones(subspcdim);
                ind=nonlin( sub2ind(size(subspc), subs), size(subspc) )*dt;
                
                if size(ind) ~= size(zz)
                    fprintf('\n\n\t\t%f\n\n',ind)
                end
                
                rrnxt(zz) = ind;
			end
		rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
		if (tspnext >= rrcum(end)) % No spike in this window
			jbin = iinxt(end)+1;
			rprev = rrcum(end);
		else   % Spike!
			ispk = iinxt(min(length(iinxt), min(find(rrcum>=tspnext)))); % time bin where spike occurred
			nsp = nsp+1;
	
			tsp(nsp,1) = ispk*dt; % spike time
			mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
			iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
			if ~isempty(iiPostSpk)
				Vmem{1}(iiPostSpk) = Vmem{1}(iiPostSpk)+ihhi(1:mxi-ispk);
			end
			tspnext = exprnd(spkrate);  % draw next spike time
			rprev = 0; % reset integrated intensity
			jbin = ispk+1;  % Move to next bin
			% --  Update # of samples per iter ---
			muISI = jbin/nsp;
			nbinsPerEval = max(minBins, round(1.5*muISI)); 
		end
	end
	tsp = tsp(find(tsp>tdim,1):nsp); % prune extra zeros
	model.spikes(round(tsp))=1;
    clear tsp;	

	[success, stats] = spikeCheck(mouse.spikes, model.spikes);
	sum1=0;
	for kk=1:length(model.spikes)
		if model.spikes(kk)==mouse.spikes(kk)
			sum1=sum1+1;
		end
	end
	perc = 100*sum1/length(model.spikes);
	meanperc=meanperc+perc;
	meansuccess=meansuccess+success;
	fitspikes{aa}=model.spikes;
	fprintf('\nData spikes: %i \nRecreated Model Spikes: %i', sum(mouse.spikes), sum(model.spikes));
	fprintf('\nPercent absolute overlap: %3.2f\n',perc);
	fprintf(['\nSucesses: %u spikes within ', num2str(1),' time bins \n\tpercent of successful data spikes: %3.2f percent \n\tpercent of successful recreated model spikes: %3.2f percent \n'],success, stats.percentTrue, stats.percentTest);
end
%}
meanperc=0;
meansuccess=0;
fitspikes=0;
%record the full conditional intensity
temp = [];
for bb=1:nvec-hasIH
    stp=[zeros(tdim-1,1); conv(stim, filt{bb}(end:-1:1),'valid')'];
    temp=[temp, stp];
end
if hasIH > 0
    stp=[zeros(length(filt{nvec})-1,1); conv(spikes, filt{nvec}(end:-1:1),'valid')' ];
    temp=[temp, stp];
end
t1=[max(1, min(BINS, round((temp-mini)/step)))];
for aa=1:length(t1)
    part1='rf(';
    part2=[];
    for bb=1:length(t1(aa,:))
        part2=strcat(part2, 't1(',num2str(aa),',',num2str(bb),')');
        if bb ~= length(t1(aa,:))
            part2=strcat(part2,',');
        end
    end
    merp=strcat(part1,part2,')');    
    condint(aa) = eval(merp); % Cond Intensity
end
%condint = rf( ind2sub(size(rf),t2) ); % Cond Intensity

meanperc=meanperc/nrep;
meansuccess=meansuccess/nrep;
