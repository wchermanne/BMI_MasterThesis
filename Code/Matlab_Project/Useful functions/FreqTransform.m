function [mySpectralData]=FreqTransform(myData,band)
%% Informations
% This functions performs a frequency analysis of the data
%
% INPUTS 
%
% myData is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% band is the string containing information about the frequency band being
% processed
%
% OUTPUTS
%
% mySpectralData is a FieldTrip Structure containing frequency information
% of myData
%
% This function also generates a figure showing the power spectrum of the
% data
%% Frequency analysis 
cfg = [];
cfg.method = 'mtmfft'; 
cfg.output = 'pow'; 
cfg.foi = [1:60];
foi=cfg.foi;
cfg.taper      = 'dpss';
cfg.tapsmofrq=2;
mySpectralData = ft_freqanalysis(cfg, myData);


%% Figures
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
