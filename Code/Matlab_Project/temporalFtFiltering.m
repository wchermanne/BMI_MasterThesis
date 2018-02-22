function [myFilteredData]=temporalFtFiltering(myFtData,LFreq,HFreq,NotchFilter,bpftype,wintype,order)
close all;
%% Preprocessing with filtering
cfg = [];
cfg.channel='all';
cfg.lpfilter      = 'yes' 
cfg.hpfilter      = 'yes' 
cfg.dftfilter     = 'yes' 
cfg.lpfreq        = HFreq;
cfg.hpfreq        = LFreq;
if(NotchFilter==1)
    cfg.dftfilter     = 'yes'
else
    cfg.dftfilter     = 'no'
end

cfg.bpfilttype=bpftype;
cfg.bpfiltwintype=wintype;
cfg.bpfiltord=order;
myFilteredData = ft_preprocessing(cfg, myFtData)

time=myFtData.time{1};
channels=myFilteredData.trial{1};
%% Figures
figure;
ax1=subplot(3,1,1)
plot(time,channels(1,:))
xlabel('time [s]')
title(['C3 x time filtered between ' num2str(LFreq) ' Hz and ' num2str(HFreq) ' Hz'])
ax2=subplot(3,1,2)
plot(time,channels(2,:))
xlabel('time [s]')
title(['C4 x time filtered between ' num2str(LFreq) ' Hz and ' num2str(HFreq) ' Hz'])
ax3=subplot(3,1,3)
plot(time,channels(3,:))
xlabel('time [s]')
title(['Cz x time filtered between ' num2str(LFreq) ' Hz and ' num2str(HFreq) ' Hz'])
end