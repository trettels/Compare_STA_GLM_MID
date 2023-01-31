ssize=get(0,'ScreenSize');
ttk = [-nkt+1:0]*reduce/10;
ff1=figure('OuterPosition',[0,0,ssize(3),ssize(4)]);
plot(ttk,gg.k/norm(gg.k),ttk,mid1.filt{1}/norm(mid1.filt{1}),ttk,mid2.filt{1}/norm(mid2.filt{1}), ttk, sta/norm(sta))
legend('GLM','MID 1d','MID + Post Spike current','STA','Location','NorthWest');
xlabel('Time preceding spike (ms)');
saveas(ff1,strcat('./figs/filt',num2str(aa),'-',num2str(bb),'.png'))

ff2=figure('OuterPosition',[0,0,ssize(3),ssize(4)]);
if (gg.iht(end) < nkt)
    newiht=gg.iht.*reduce/10;
else
    endex=find(round(gg.iht)>=nkt,1);
    newiht=gg.iht(1:endex);
end
termin=round(newiht(end));
glmih=gg.ihbas*gg.ih;
glmih=glmih/norm(glmih);
newiht=newiht*reduce/10;
mid2.ih=resample(mid2.filt{2}(1:termin),length(newiht),termin, 2);
mid2.ih=mid2.ih/norm(mid2.ih);
plot(newiht,glmih(1:endex),newiht,mid2.ih)
legend('glm','mid + post-spike')
title('Post-spike Current')
xlabel('Time after spike (ms)')
saveas(ff2,strcat('./figs/ih',num2str(aa),'-',num2str(bb),'.png'))

ff3=figure('OuterPosition',[0,0,ssize(3),ssize(4)]);
hist(diff(tsp*reduce/10),100)
title('ISI')
xlabel('Interspike Interval (ms)');
saveas(ff3,strcat('./figs/isi',num2str(aa),'-',num2str(bb),'.png'))

ff4=figure('OuterPosition',[0,0,ssize(3),ssize(4)]);
batman=diff(tsp)*reduce/10;
hist(batman(find(batman<=nkt*reduce/10)),100)
title('ISI')
xlabel('Interspike Interval (ms)');
saveas(ff4,strcat('./figs/isi',num2str(aa),'-',num2str(bb),'_zoom.png'))