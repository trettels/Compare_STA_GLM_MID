%%Generator v3: Powell's direction set method
function[newvec]=genParPowell(vec, stim,spikes,nfilt, BINS, nvec, tol)
global HCURRENT;
nvec=round(nvec);
if ~isscalar(nfilt)
    nfilt=nfilt(end);
end

fprintf('=');
if(nvec > 1 && HCURRENT > 0)
	rr=1;
else
	rr=0;
end
for kk=1:nvec-rr
	vec{kk}=-vec{kk}/norm(vec{kk});
end

newvec=cell(1,nvec);

% Call powell's method
[pt]=testPowell(vec, nfilt, nvec, stim, spikes, BINS, tol);

%  Reassign to proper ouput
for p=1:nvec
	newvec{p}=pt(:,p)';
end
for(p=1:nvec-rr)
    newvec{p}=newvec{p}/norm(newvec{p});
end
