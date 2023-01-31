function [testPSTH, inputPSTH, fig] = rateCodingCheck(test, inputPSTH, stim, bins, nrep)
%rateCodingCheck.m
setpaths
global BINSIZE;
BINSIZE=bins;
testPSTH=zeros(1,length(stim));
global HCURRENT;
hasIH = HCURRENT;
t_nvec=length(test.filt);

if isfield(test,'rf') && ~isfield(test,'rf1')
    test.rf1=test.rf;
end
   
[ts] = getTrain(length(test.filt{1}),stim, test.filt, test.rf1, length(test.filt), 1, nrep);

for aa=1:nrep
    testPSTH=testPSTH + ts{aa};
end

testPSTH=testPSTH./nrep;

ssize=get(0,'ScreenSize');
fig=figure('OuterPosition',[0,0,ssize(3),ssize(4)]); clf;
kl=length(test.filt{1});
plot(1+kl:1000+kl, inputPSTH(1+kl:1000+kl),'-rx', 1+kl:1000+kl, testPSTH(1+kl:1000+kl),'-bo')
legend('Model', 'Test (Derrived model)');


