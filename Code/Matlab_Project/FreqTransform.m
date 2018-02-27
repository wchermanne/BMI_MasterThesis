function [mySpectralData]=FreqTransform(data,band)
%% Frequency analysis 
cfg = [];
cfg.method = 'mtmfft'; 
cfg.output = 'pow'; 
cfg.foi = [1:60];
foi=cfg.foi;
cfg.taper      = 'dpss';
cfg.tapsmofrq=2;
mySpectralData = ft_freqanalysis(cfg, data);
freq=mySpectralData.freq;
power_spectrum=mySpectralData.powspctrm;

titleC3='Power spectrum of C3 for ';
myTitleC3=[titleC3 band];
titleC4='Power spectrum of C3 for ';
myTitleC4=[titleC4 band];
titleCz='Power spectrum of C3 for ';
myTitleCz=[titleCz band];

figure;
ax1 = subplot(3,1,1);
plot(foi,power_spectrum(1,:))
xlabel('time [s]')
ylabel('power')
title(myTitleC3)

ax2 = subplot(3,1,2);
plot(foi,power_spectrum(2,:))
xlabel('time [s]')
ylabel('power')
title(myTitleC4)

ax3 = subplot(3,1,3);
plot(foi,power_spectrum(3,:))
xlabel('time [s]')
ylabel('power')
title(myTitleCz)
end
