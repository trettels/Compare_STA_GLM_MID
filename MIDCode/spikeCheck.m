%  proj: Vector of stimulus projections onto the filter, run through the
%  nonlinearity
%  
%  thresh: spike threshold

function [success, stats] = spikeCheck(truespikes, testspikes)

%Initializations
success=0;
lastspiketime=-1;
lastpostime=-1;
winsize = 1;
for i=1:length(truespikes)
    
    % Check to see if we should discard the last lnp spike time
    if i-lastpostime > winsize && lastpostime ~=-1
        lastpostime=-1;
    end
    
    %check to see if we should discard the last model spike time
    if i-lastspiketime > winsize && lastspiketime ~=-1
        lastspiketime=-1;
    end
    
    % Are there two perfectly-overlapping spikes?
    if truespikes(i)==1 && testspikes(i)==1
        success = success+1;
        lastspiketime=-1;
        lastpostime=-1;
    % Is there a model spike to store?
    elseif truespikes(i)==1 && lastpostime ==-1 && testspikes(i)==0
        lastspiketime=i;
    % Is there an poisson spike to store?
    elseif testspikes(i)==1 && lastspiketime == -1 && truespikes(i)==0
        lastpostime=i;
    % Is there a prior poisson spike within the window of the model spike?
    elseif truespikes(i)==1 && lastpostime ~=-1 && (i-lastpostime <= winsize || testspikes(i)==1)
        success = success+1;
        lastspiketime=-1;
        lastpostime=-1;
    % Is there a model spike within the window of the poission spike?
    elseif testspikes(i)==1 && lastspiketime ~= -1 && i-lastspiketime<=winsize
        success=success+1;
        lastpostime =-1;
        lastspiketime=-1;
    end
end

stats.percentTrue = 100*success/sum(truespikes);
stats.percentTest = 100*success/sum(testspikes);
