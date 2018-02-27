function [comp] = ICA(myData)
% Laplacian performs a spatial filtering following the Independent Component Analysis method

% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
%
% In the param structure, the sample frequency is saved
% 
% The ICA function returns the filtered EEG signals
channels=myData.data(1).trial{1};
time=myData.data(1).time{1};
cfg=[];
cfg.method='runica';
comp=ft_componentanalysis(cfg,myData.data(1));

channels=comp.trial{1};
figure;
ax1=subplot(3,1,1);
plot(time,channels(1,:))
xlabel('time [s]')
title('Comp channel 1')
ax2=subplot(3,1,2);
plot(time,channels(2,:))
xlabel('time [s]')
title('Comp channel 2')

ax3=subplot(3,1,3);
plot(time,channels(3,:))
xlabel('time [s]')
title('Comp channel 3')


end
