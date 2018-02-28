function [myFilteredData]=temporalFtFiltering(myFtData,LFreq,HFreq,NotchFilter,bpftype,wintype,order,band)
%% Informations
% This functions preprocesses the data
%
% INPUTS 
%
% myFtData is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% LFreq is the low pass frequency of the BP filter, in Hz
%
% HFreq is the high pass frequency of the BP filter, in Hz
%
% NotchFilter is 1 if a NotchFilter needs to be applied, 0 if not
%
% bpftype is the type of BP filter used. 'firws','but' or 'fir'
%
% wintype is the type of window used for firws filters
%
% order is the order of the bp filter 
%
% band is the string containing information about the frequency band being
% processed
%
% OUTPUTS
%
% data_filtered is the FieldTrip Structure containing the filtered data


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

titleC3='C3 x time temporal filtered for '
myTitleC3=[titleC3 band];

titleC4='C4 x time temporal filtered for '
myTitleC4=[titleC4 band];

titleCz='Cz x time temporal filtered for '
myTitleCz=[titleCz band];


%% Figures
figure;
ax1=subplot(3,1,1)
plot(time,channels(1,:))
xlabel('time [s]')
title(myTitleC3)
ax2=subplot(3,1,2)
plot(time,channels(2,:))
xlabel('time [s]')
title(myTitleC4)
ax3=subplot(3,1,3)
plot(time,channels(3,:))
xlabel('time [s]')
title(myTitleCz)
end