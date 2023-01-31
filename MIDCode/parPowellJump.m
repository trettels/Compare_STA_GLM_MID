function [newVec] = parPowellJump(testVec, T, stim, spikes, nfilt_p, BINS, nvec, tol)
nvec=round(nvec);
global HCURRENT;
%HCURRENT=0;
%Ensure the integers are actually integers...
nfilt_p=round(nfilt_p);
if ~isscalar(nfilt_p)
    nfilt=prod(nfilt_p);
else
    nfilt=nfilt_p;
end

BINS=round(BINS);
vec=zeros(nfilt,nvec);
if iscell(testVec)
    for k=1:nvec
        vec(:,k)=testVec{k};
    end
end

% Random offset
ranseed=.1-.2*rand(nfilt,nvec);
if nvec > 1
    [ranseed,~,~]=svd(ranseed,0);
else
    ranseed=ranseed/norm(ranseed);
end  
vec=vec+T*ranseed;

%{
%Normalize vectors
for kk=1:nvec
    vec(:,kk)=vec(:,kk)/norm(vec(:,kk));
end
%}
%Zero-out Ih, if present
%if nvec > 1 && HCURRENT==1
%    vec(:,end) = ones(nfilt,1);
%end

% Go through one minimization 
%nVec=testPowell(vec, nfilt, nvec, stim, spikes, BINS, tol);

% Put everything in the proper form to return
for kk=1:nvec
    newVec{kk}=vec(:,kk)';
end
