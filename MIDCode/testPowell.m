function [pt]=testPowell(vec, nfilt, nvec, stim, spikes, BINS, TOL)
% Implements Powell's method from Numerical Recipies (Section 10.7, 3rd Ed)
% Still needs to be re-adjusted for my data structures... send in vectors
% as matrices or cell-arrays?
global HCURRENT;
if HCURRENT > 0
    ihcurr=1;
else
    ihcurr=0;
end
[m,n]=size(nfilt);
if m~=n
    error('powell:argCheck','\n"nfilt" should be a scalar value, not a matrix\n');
end
% Make sure the vector is in the proper format:
if ~iscell(vec)
    [m,n]=size(vec);
    if m<n
        pp=vec';
        fprintf('\nInput vector improperly aligned, taking the transpose...\n');
    else
        pp=vec;
    end
else
    pp=zeros(nfilt,nvec);
    for k=1:nvec  %-ihcurr
        pp(:,k)=vec{k}'/norm(vec{k});
    end
end

% Set up the constants
MAXIT=200;
%TOL=1e-3; %1e-10;
%fprintf('\nTolerance: %e\n', TOL);
LINTOL=1e-10;
LINITER=200;
% Set up the variables
p=pp;
n=nfilt*nvec; % The "vector" length / iteration length
ximat=eye(n);  %The initial direction set
fret=findNegInfoMEX(stim,spikes,p, nfilt, BINS,nvec);
pt=p;

xi=cell(1,n);
for iter=1:MAXIT
    %fprintf('.');
    %Eval the current direction set
    fp=fret;
    oldInfo = repmat(fp,1,n);
    delta=zeros(1,n);
    del=0.0;
    xii=zeros(nfilt,nvec);
    parfor i=1:n
        xi{i}=ximat(:,i);
        % reshape xi into a pair of vectors
        xi{i}=(reshape(xi{i},nfilt,nvec));

        %fptt=fret;
        %Update p
        V=cell(1,nvec);
        X=cell(1,nvec);
        for k=1:nvec
            V{k} = (p(:,k))'; %+updater(:,k)'; % Use a filter set which has been updated with the local displacement vectors
            X{k}=xi{i}(:,k)';
        end
        %fprintf('starting line minimization\n');
        lambda=linearMinimizer(V,X,LINTOL,LINITER,stim,spikes,nfilt,BINS, nvec, ihcurr);
        for k=1:nvec
            xi{i}(:,k)=xi{i}(:,k)*lambda(k);
        end
        
%       updater=updater+xi{i};
        %p=p+xi{i};
        %{
        labBarrier;

        xii=xii+gplus(updater); %Update p across labs

        g{i}=g{i}+updater; % Local update to the test point, to assist in convergance
        labBarrier; 
        %} 
        info(1,i) = findNegInfoMEX(stim,spikes,p+xi{i}, nfilt, BINS,nvec, ihcurr);   %p+xi{i} 
    end

    %fprintf('finished minimizing iteration');
    upd=zeros(nfilt,nvec);

    parfor zz=1:n
        upd(zz)=xi{zz}(zz);
    end
    pz=p;
    p=p+upd;

%{
    % Normalize the stimulus filters
    ptem=p(:,1:nvec-HCURRENT);
    ih=p(:,end);
    if nvec-HCURRENT > 1
        fprintf('^');
        [ptem,~] = svd(ptem,0);
    else
        ptem=ptem/norm(ptem);
    end
    if HCURRENT > 0
        fprintf('!');
        p=[ptem, ih];
    else
        p=ptem;
    end
%}      
    delta = oldInfo-info;
    ibig=find(abs(delta)==max(abs(delta)),1);  % ensure find only returns one index (the first one)
    del=delta(ibig);
    % Normalize Ih if it is too big
    for rr=1:nvec  %-ihcurr
    	    p(:,rr) = p(:,rr)/norm(p(:,rr)); %zeros(size(p(:,2)));
        	%fprintf('#');
    end

    fret=findNegInfoMEX(stim,spikes,p, nfilt, BINS,nvec);
    if 2*(fp-fret)<=TOL*(abs(fp)+abs(fret))+eps %|| fp < fret
        if(fret <= fp) % || iter==1)
            %fprintf('\nold info: %f, new info: %f\n',-fp,-fret);
            pt=p;
        elseif(iter==1)
            pt=pt+TOL*(rand(size(pt))*2-1);
        end
        break
    elseif iter==MAXIT
        fprintf('Exceeding max iterations in Powell\n');
        pt=p;
    end
    
    ptt=2*p-pt;
    xi_new=p-pt;
    pt=p;
    fptt=findNegInfoMEX(stim,spikes,ptt, nfilt, BINS,nvec);
    V=cell(1,nvec);
    X=cell(1,nvec);
    if fptt < fp
        t=2*(fp-2*fret+fptt)*(fp-fret-del)^2-del*(fp-fptt)^2;
        if t<0.0
            %Update p
            for k=1:nvec
                V{k}=p(:,k)';
                X{k}=xi_new(:,k)';
            end
            lambda_2=linearMinimizer(V,X,LINTOL,LINITER,stim,spikes,nfilt,BINS, nvec, ihcurr);
            for k=1:nvec
                p(:,k)=p(:,k)+xi_new(:,k)*lambda_2(k);
            end
%            [p,~,~]=svd(p,0);
			for rr=1:nvec  %-ihcurr
				p(:,rr) = p(:,rr)/norm(p(:,rr)); %zeros(size(p(:,2)));
				%fprintf('#');
			
			end
       %{     
            % Normalize and orthoganalize (if relevant)
            ptem=p(:,1:nvec-HCURRENT);
            ih=p(:,end);
            if nvec-HCURRENT > 1
                [ptem,~] = svd(ptem,0);
            else
                ptem=ptem/norm(ptem);
            end
            if HCURRENT > 0
                p=[ptem, ih];
            else
                p=ptem;
            end
       %}     
            fret=findNegInfoMEX(stim,spikes,p, nfilt, BINS,nvec);
%{
            % Determine where the direction goes
            if ibig <=n/2
                inser = n/2;
            else
                inser=n;
            end
%}            
            inser=n;
            ximat(:,ibig)=ximat(:,inser);
            xi_new=reshape(xi_new,n,1);
            ximat(:,n)=xi_new;
        end
    end

    %p(:,2)=p(:,2)/norm(p(:,2));       

end

%fprintf('\nExiting Powell\n');
