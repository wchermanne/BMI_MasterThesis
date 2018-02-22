clear all; close all;

%% Load the files and data
load('parameters.mat'); % Parameters of eeg_generator
load('rawData.mat'); % Data generated by eeg_generator
rawData=[C3;C4;Cz]


%% Create the FieldTrip structure
fsample=parameters.fsample;
t_end=parameters.duration;

myFtData=[];
myFtData.fsample=fsample;
myFtData.trial={rawData};
myFtData.label={'C3' 'C4' 'CZ'};
myFtData.sampleinfo=[1 length(C3)];
myFtData.time=parameters.time;
myFtData.cfg=[];


%% Preprocessing without filtering
cfg = [];
cfg.channel='all';
[pre_processed_data] = ft_preprocessing(cfg, myFtData);


%% Show data
%[cfg] = ft_databrowser(cfg, pre_processed_data)


%% Frequency analysis 
cfg = [];
cfg.method = 'mtmfft'; 
cfg.output = 'pow'; 
cfg.foi = [1:60];
cfg.taper      = 'dpss';
cfg.tapsmofrq=2;
freqdata = ft_freqanalysis(cfg, pre_processed_data);
freq=freqdata.freq;
power_spectrum=freqdata.powspctrm;


%% Preprocessing with filtering
cfg = [];
cfg.channel='all';
cfg.lpfilter      = 'yes' 
cfg.hpfilter      = 'yes' 
cfg.dftfilter     = 'yes' 
cfg.lpfreq        = 12;
cfg.hpfreq        = 8;
lpfreq=cfg.lpfreq;
hpfreq=cfg.hpfreq;

[pre_processed_data_filter] = ft_preprocessing(cfg, mydata)

channels_filtered= pre_processed_data_filter.trial{1};
channels_filtered_squared=channels_filtered.^2;

%% Show data
%[cfg] = ft_databrowser(cfg, pre_processed_data)


%% Frequency analysis 
cfg = [];
cfg.method = 'mtmfft'; 
cfg.output = 'pow'; 
cfg.foi = [1:60];
foi=cfg.foi;
cfg.taper      = 'dpss';
cfg.tapsmofrq=2;
freqdata_filter = ft_freqanalysis(cfg, pre_processed_data_filter);
freq=freqdata_filter.freq;
power_spectrum_filtered=freqdata_filter.powspctrm;
cfg=[];

%% Figures and plots

%% Signals at input
figure;
ax1 = subplot(3,1,1);
plot(time,channels(1,:))
xlabel('time [s]')
title('C3 x time at the input')
ax2=subplot(3,1,2)
plot(time,channels(2,:))
xlabel('time [s]')
title('C4 x time at the input')
ax3=subplot(3,1,3)
plot(time,channels(3,:))
xlabel('time [s]')
title('Cz x time at the input')

%% Signals after preprocessing

figure;
ax1=subplot(3,1,1)
plot(time,channels_filtered(1,:))
xlabel('time [s]')
title(['C3 x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])
ax2=subplot(3,1,2)
plot(time,channels_filtered(2,:))
xlabel('time [s]')
title(['C4 x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])
ax3=subplot(3,1,3)
plot(time,channels_filtered(3,:))
xlabel('time [s]')
title(['Cz x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])


%% Signals after squaring

figure;
ax1 = subplot(3,1,1);
plot(time,channels_filtered_squared(1,:))
xlabel('time [s]')
ylabel('power')
title(['C3 x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz and squared'] )

ax2 = subplot(3,1,2);
plot(time,channels_filtered_squared(2,:))
xlabel('time [s]')
ylabel('power')
title(['C4 x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz and squared'] )

ax3 = subplot(3,1,3);
plot(time,channels_filtered_squared(3,:))
xlabel('time [s]')
ylabel('power')
title(['Cz x time filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz and squared'] )

%% Signals in frequency domain

figure;
cfg=[];
ft_singleplotER(cfg,freqdata)
xlabel('Freq (Hz)')
ylabel('Power')
title('Power spectrum of the non-filtered data (mean of channels)')

figure;
ax1 = subplot(3,1,1);
plot(foi,power_spectrum(1,:))
xlabel('time [s]')
ylabel('power')
title('Power spectrum of C3 non-filtered')

ax2 = subplot(3,1,2);
plot(foi,power_spectrum(2,:))
xlabel('time [s]')
ylabel('power')
title('Power spectrum of C4 non-filtered')

ax3 = subplot(3,1,3);
plot(foi,power_spectrum(3,:))
xlabel('time [s]')
ylabel('power')
title('Power spectrum of Cz non-filtered')

figure;
cfg=[];
ft_singleplotER(cfg,freqdata_filter)
xlabel('Freq (Hz)')
ylabel('Power')
title(['Power spectrum of the data filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'] )

figure;
ax1 = subplot(3,1,1);
plot(foi,power_spectrum_filtered(1,:))
xlabel('time [s]')
ylabel('power')
title(['Power spectrum of C3 filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])

ax2 = subplot(3,1,2);
plot(foi,power_spectrum_filtered(2,:))
xlabel('time [s]')
ylabel('power')
title(['Power spectrum of C4 filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])

ax3 = subplot(3,1,3);
plot(foi,power_spectrum_filtered(3,:))
xlabel('time [s]')
ylabel('power')
title(['Power spectrum of Cz filtered between ' num2str(hpfreq) ' Hz and ' num2str(lpfreq) ' Hz'])
