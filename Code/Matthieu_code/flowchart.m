%%% Flowchart with different operations required for conditionning the BCI
%%% signals

%% Signal generation for test (continuous sensorimotor rythm ) using 1st Order MPA
close all;
Fs = 250; % sampling rate
freq = linspace(0,35,36);
freq = [freq 50 60]';

amp = [0.5 4/8 7/8 4/8 2.5 0.3 0.15 0.3 3 5 6 8 5 4 2 1.5 1.2 1 0.9 0.6 0.4 0.2 0.1 0.1 0.05 0.03 0.02 0.01 0.008 0.004 0.002 0.001 0.0008 0.0005 0.0003 0.0001 20 20]*10;
phi = pi/2;
t_end = 15;
t = 0:1/Fs:t_end
sin_mat = zeros(1,length(t));
logbook_gamma = zeros(1,length(t));
logbook_amp =  zeros(1,length(t));
desync_duration = 0.4 + 0.1*rand;
time_movement = 9;
movement_duration = 1;

%%% Taking eye blink effect into account
%eye_blink_occurence = floor(2*rand);
eye_blink_occurence = 1;
time_eye_blink = ones(1,eye_blink_occurence);
eye_blink_rising_duration = 0.075 + 0.01*rand;
eye_blink_steady_duration = ones(1,eye_blink_occurence);
for k =1:1:eye_blink_occurence
    time_eye_blink(k) = 1+(t_end./2)*k*rand;
    eye_blink_steady_duration(k) =0.2 + 0.2*rand;
end



for l = 1:1:length(t)
    amp = amp.*(amp>0); %% for making amp always positive
    %%% Taking Event-related desynchronization and resync for lower mu
    %%% rhythm. That happens for any motor task 
    if (t(l) <= time_movement && t(l) >= time_movement-2*desync_duration)
        amp(8:10) = 0.95.*amp(8:10);
    end
    if (t(l) <= time_movement && t(l) >= time_movement-desync_duration)
        amp(8:10) = 1.*amp(8:10)./0.95;
    end
    %%% Taking Beta rhythm 
    if (t(l) >= time_movement && t(l) <= time_movement+movement_duration)
        amp(18:26) = 1.01.*amp(18:26);
    end
    if (t(l) >= time_movement+movement_duration && t(l) >= time_movement-2*movement_duration)
        amp(18:26) = amp(18:26)./(1.01);
    end
    
    
    %%% Taking Eye Blink
    if(isempty(time_eye_blink) ==1)
        
    else
        for k =1:1:eye_blink_occurence
            if(t(l) >= time_eye_blink(k) && t(l) <=time_eye_blink(k)+eye_blink_rising_duration)
                amp(1:2) = 1.3.*amp(1:2);
            end
            if(t(l) >= time_eye_blink(k) + eye_blink_rising_duration +eye_blink_steady_duration && t(l) <=time_eye_blink(k)+2*eye_blink_rising_duration+eye_blink_steady_duration)
                amp(1:2) = 1.*amp(1:2)./(1.3);
            end
        end
        
    end
    logbook_amp(l) = amp(10);
    sin_mat_temp = 0;
    for k = 1:1:length(freq)
        sin_mat_temp = sin_mat_temp + amp(k)*sin(2*pi*freq(k)*t(l)+phi);
    end
    sin_mat(l) = sin_mat_temp;
    gamma = 0.99 + 0.02*rand;
    logbook_gamma(l) = gamma;
    amp = gamma*amp; + 0.5*rand(1,length(amp)) - 0.25;
    amp(end-1:end) = [200 200];
end
figure; plot(logbook_amp)
%%% Add noise contribution
Y = awgn(sin_mat,10,'measured');

figure;
subplot(2,1,1)
plot(t,Y)

dft_size = 512;
fourier_sig = fftshift(fft(Y,dft_size));
k = 0:1:dft_size-1;
f_axis = Fs*k/dft_size -Fs/2;
fourier_sig = mag2db(abs(fourier_sig));
subplot(2,1,2)
plot(f_axis,fourier_sig);
xlabel('frequency [Hz]');
ylabel('Magnitude [dB]');
grid on

nyq_freq = Fs./2; %% Half the sampling rate; nyquist frequency
wlen = 512;
hop = 256;
[stft_tot, freq_vec, time_vec] = stft(Y, wlen, hop, dft_size, Fs);
figure;
surf(time_vec,freq_vec',mag2db(abs(stft_tot)));

figure;
spectrogram(Y,512,256,dft_size,Fs,'power');
%% Frequency Filtering
%%% One might use bandpass filtering from 3Hz (larger than EOG and eye blink artifacts)
%%% up to 30-40 Hz (lower than power line interference @50 Hz)
low_cutoff = 8; % 3 Hz
high_cutoff = 13; % 35Hz

%%% 1ST : FIR Filter --> pay attention to the delay = N/2
b_fir = fir1(192*2,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
figure;
freqz(b_fir,1,512);

%%% 2ND : IIR Filter. 2 possibilities : Either Butterworth filter
%%% (maximally flat filter) or Chebychev 2 (flat in the pass band and
%%% distortion in the stop band)
[b_iir,a_iir]=cheby2(8,60,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]); % Bandpass digital filter design (12th order, 60dB attenuation in the stopband)
%[b_iir,a_iir] = zp2tf(z_iir,p_iir,k_iir);
%iir_filt = tf(b_iir,a_iir,1/Fs);
%pzmap(iir_filt);
%h = fvtool(b_iir,a_iir);

%%% One might alos compare it to an elliptic filter or cheby1


%%% Results obtained with FIR filter
Y_filtered = filter(b_fir,1,Y)
figure
plot(t,Y_filtered)
figure
dft_size = 512;
fourier_sig = fftshift(fft(Y_filtered,dft_size));
k = 0:1:dft_size-1;
f_axis = Fs*k/dft_size -Fs/2;
fourier_sig = mag2db(abs(fourier_sig));
plot(f_axis,fourier_sig);
xlabel('frequency [Hz]');
ylabel('Magnitude [dB]');
grid on

%%% STFT after filtering the signal :
wlen =512;
hop = 256;
[stft_tot, freq_vec, time_vec] = stft(Y_filtered, wlen, hop, dft_size, Fs);
figure;
surf(time_vec,freq_vec',mag2db(abs(stft_tot)));
figure; plot(time_vec,abs(stft_tot(26,:)),time_vec,abs(stft_tot(25,:)),time_vec,abs(stft_tot(24,:)))
%% Spatial Filtering

%% Windowing and Framing
% initialize the signal time segment index
indx = 0;
winlen = 100; %% in ms
winoverlap = 50; %% in ms
%wlen = ((winlen./1000)*Fs); %% window length in samples
wlen = 24;
%hop = floor((winoverlap./1000)*Fs); %% window overlap in samples
hop = 12;
win = hamming(wlen, 'periodic');

%%% Need to compute coln /!\ what about resamping ? simplest idea is to get
%%% rid of the last so many samples so that you obtain an integrer number
frame_numb = floor(length(Y_filtered)./hop);
frame_numb_removing= length(Y_filtered)-(floor(length(Y_filtered)./hop))*hop;
Y_filtered = Y_filtered(1:end-frame_numb_removing);
t = t(1:end-frame_numb_removing);
coln = frame_numb;
Yw_matrix = zeros(wlen,coln./2);
for k = 1:1:coln-1
    % windowing
    Yw = (Y_filtered(indx+1:indx+wlen).').*win;
    % Matrix with time segements
    Yw_matrix(:,k) = Yw;
    % update the index
    indx = indx + hop;
end

%% Feature extraction - Time domain
%%% 1st method : peak detection
peak = zeros(coln-1,1);
for k = 1:1:coln-1
    peak(k) = max(Yw_matrix(:,k));
end
figure;
subplot(2,1,1)
plot(t,Y_filtered);
subplot(2,1,2);
plot(1:1:coln-1,peak);
title('Peak Picking Method')

%%% 2nd method : Integration. This consists of squaring the signal and
%%% integrate it. (Integration is made using Darboux sum or trapezoidal integration. Using higher order methods is not relevant here since
%%% we aren't looking for the exact value. Notice Darboux sum is just a
%%% scalar product
integration = zeros(coln-1,1);
for k = 1:1:coln-1
    integration(k) = ((Yw_matrix(:,k)).')*((1/Fs)*ones(wlen,1));
end
figure;
subplot(2,1,1)
plot(t,Y_filtered);
subplot(2,1,2);
plot(1:1:coln-1,integration);
title('Inegration Method')

%% Feature extraction - Frequency Domain
%%% 1st method : Band power --> Rather equivalent to time integration
integration = zeros(coln-1,1);
for k = 1:1:coln-1
    integration(k) = ((Yw_matrix(:,k).^2).')*((1/Fs)*ones(wlen,1));
end
figure;
subplot(2,1,1)
plot(t,Y_filtered);
subplot(2,1,2);
plot(1:1:coln-1,integration);
title('Band Power Method')

%%% 2nd method : This is equivelent to perform a STFT (one could explain the algorithm)
% rown = ceil((1+dft_size)/2);
% stft_mat = zeros(rown,coln-1);%frequency across the colums and time across the rows
% for k = 1:1:coln-1
%     X = fft(Yw_matrix(:,k),dft_size);
%     stft_mat(:,k) = X(1:rown);
% end
% t_stft = (wlen/2:hop:wlen/2+(coln-1)*hop)/Fs;
% f_stft = (0:rown-1)*Fs/dft_size;
% figure;
% surf(t_stft,f_stft',mag2db(abs(stft_mat)));

%%% 3nd method : Cepstra. This is a Laurents serie of the spectrum. We know
%%% that the spectrum is periodic --> can be expressed as a truncated
%%% cosine serie


%%% 4nd method : Linear predictive coding, also called Autoregressive
%%% model. The current sample might be predicted aq a linear combination of
%%% the previous samples + a prediction error --> Levinson Durbin algorithm
%%% = O(n^2) instead of O(n^3) with covariance method
p = 10; % order of the AR model
r=zeros(p+1,size(Yw_matrix,2));
% Yw_matrix = (wlen-by-frame_number-1) framing matrix, one column per frame of
% length wlen
% p = number of lags
% r = (p+1)-by-frame_number of unnormalized autocorrelations, r(0,:) on top, r(p,:)
%      at the bottom
for k=1:p+1
    r(k,:)=sum(Yw_matrix(k:end,:).*Yw_matrix(1:end-k+1,:),1);
end

% function [a,k,Eres]=cor2lpc(r)
% a = predicition polynomial coefficients, including the leading 1 (a0)
% Beware for sign: use polynomial coefficients as-is: freqz(1,a(:,frame))
% is transfer function
% k = reflection coefficients
% Eres: variance of the residual

[N,T]=size(r);
a=zeros(N,T);a(1,:)=1;
k=zeros(N-1,T);

for m=1:N-1,
    k(m,:)=sum(r(m+1:-1:2,:).*a(1:m,:),1) ./ sum(r.*a,1);
    a(2:m+1,:)=a(2:m+1,:)-k(m*ones(m,1),:).*a(m:-1:1,:);
end
Eres=sum(r.*a,1);

%% Feature Extraction - Mixed Frequency/Time methods
