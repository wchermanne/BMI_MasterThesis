feature_mat = ones(12,3);
load('LDA_with_CSP.mat')
for i = 12:1:12
    %% Loading data
    load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/Data_for_CSP/EEG_Signals_Left_Trial_' num2str(i) '.mat']);
    nyq_freq = Fs/2;
    
    %% Frequency Filtering & subband creation for SBCSP with 4 subbands
    %%% Here a filterbank is designed going from 7 to 31Hz
    %%% up to 30-40 Hz (lower than power line interference @50 Hz)
    filter_order = 96;
    
    %%% 1ST : FIR Filter for 7-> 13Hz pay attention to the delay = N/2
    low_cutoff = 7; % Hz
    high_cutoff = 13; % Hz
    b_fir_1 = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %%% 2ND : FIR Filter for 13-> 19Hz pay attention to the delay = N/2
    low_cutoff = 13; % Hz
    high_cutoff = 19; % Hz
    b_fir_2 = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %%% 3RD : FIR Filter for 19-> 25Hz pay attention to the delay = N/2
    low_cutoff = 19; % Hz
    high_cutoff = 25; % Hz
    b_fir_3 = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %%% 4TH : FIR Filter for 25-> 31Hz pay attention to the delay = N/2
    low_cutoff = 25; % Hz
    high_cutoff = 31; % Hz
    b_fir_4 = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %% Temporal Filtering
    numOfBand = 4;
    total_fir = [b_fir_1;b_fir_2;b_fir_3;b_fir_4];
    Channel.time = time;
    Channel.Fs = Fs;
    for k=1:1:numOfBand
        actual_filter= total_fir(k,:);
        C3_filtered = filter(actual_filter,1,C3); % C3_prime but let's modify it
        C4_filtered = filter(actual_filter,1,C4);
        Cz_filtered = filter(actual_filter,1,Cz);
        Channel.band(k).C3 = C3_filtered;
        Channel.band(k).C4 = C4_filtered;
        Channel.band(k).Cz = Cz_filtered;
    end
    
    
    %% Multi-Band CSP + Feature Vector generation
    
    %%% 3RD : Common spetcral patteren %%%. Notice that this is an adaptive
    %%% technique
    % Notice that we assume that X = [C3 Cz C4] (from left to right)
    for k=1:1:numOfBand
        R_ave_right = zeros(3,3);
        for trial = 1:1:10
            load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Right_Trial_' num2str(trial) '.mat']);
            R_normalized_current = spatial_cov_computation(Channel.Fs,Channel.band(k).C3,Channel.band(k).C4,Channel.band(k).Cz,Channel.time);
            R_ave_right = R_ave_right + R_normalized_current;
        end
        R_ave_right = R_ave_right./10;
        
        R_ave_left = zeros(3,3);
        for trial = 1:1:10
            load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Left_Trial_' num2str(trial) '.mat']);
            R_normalized_current = spatial_cov_computation(Channel.Fs,Channel.band(k).C3,Channel.band(k).C4,Channel.band(k).Cz,Channel.time);
            R_ave_left = R_ave_left + R_normalized_current;
        end
        R_ave_left = R_ave_left./10;
        R_tot = R_ave_right + R_ave_left;
        [U,D] = eig(R_tot);
        P = (inv(D))^(1/2)*(U.');
        Sigma_1_hat = P*R_ave_right*(P.');
        Sigma_2_hat = P*R_ave_left*(P.'); % The sum of the 2 previoulsy
        [V,Gamma] = eig(Sigma_1_hat);
        W = (P.')*V;
        %%%% Then we can select first and last vector. They have max variance
        %%%% for class 1 & 2 respectively
        W_Csp = [W(:,1) W(:,2)];
        A = (W.')*[Channel.band(k).C3; Channel.band(k).Cz; Channel.band(k).C4];
        S.band(k).Ch1 = A(1,:);
        S.band(k).Ch2 = A(2,:);
        S.time = Channel.time;
        S.Fs = Channel.Fs;
        %%%% FEATURE VECTOR GENERATION AROUND T = 10 (MOVEMENT TIME) %%%%%%
        t_move = 10;
        index_time = find(time ==t_move);
        wlen = 3*Fs;
        win = (hamming(wlen, 'periodic')).';
        S1_win = S.band(k).Ch1(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
        S2_win = S.band(k).Ch2(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
        var_tot = var(S1_win) + var(S2_win); %% Used for normalisation !
        S.feature_vectors(k).cmpt1 = log(var(S1_win)./var_tot);
        S.feature_vectors(k).cmpt2 = log(var(S2_win)./var_tot);
    end
end
close all

%%% First trial of machine learning prediction






%%


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

figure;
spectrogram(Y_filtered,512,256,dft_size,Fs,'power');
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
coln = frame_numb- frame_numb_removing;
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
