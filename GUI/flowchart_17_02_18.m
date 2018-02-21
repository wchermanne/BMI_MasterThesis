feature_mat = ones(12,4);
load('LDA_2.mat')
for i = 12:1:12
    %% Loading data
    load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/Data_for_CSP/EEG_Signals_Right_Trial_' num2str(i) '.mat']);
    nyq_freq = Fs/2;
    Y = C3;
    t = time;
    figure;
    subplot(3,1,1)
    plot(time,C3)
    subplot(3,1,2)
    plot(time,Cz)
    subplot(3,1,3)
    plot(time,C4)
    
    %% Frequency Filtering & suuband creation for
    %%% One might use bandpass filtering from 3Hz (larger than EOG and eye blink artifacts)
    %%% up to 30-40 Hz (lower than power line interference @50 Hz)
    filter_order = 96;
    
    %%% 1ST : FIR Filter for mu rhythm --> pay attention to the delay = N/2
    low_cutoff = 7; % 3 Hz
    high_cutoff = 14; % 35Hz
    b_fir_mu = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    %figure;
    %freqz(b_fir_mu,1,512);
    
    %%% 2ND : FIR Filter for beta rhythm --> pay attention to the delay = N/2
    low_cutoff = 15; % 3 Hz
    high_cutoff = 26; % 35Hz
    b_fir_beta = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    %figure
    %freqz(b_fir_beta,1,512);
    
    %%% 3RD : FIR Filter for gamma rhythm --> pay attention to the delay = N/2
    low_cutoff = 35; % 3 Hz
    high_cutoff = 45; % 35Hz
    b_fir_gamma = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    %figure
    %freqz(b_fir_gamma,1,512);
    
    %%% ALT FILTER : FIR Filter for gamma rhythm --> pay attention to the delay = N/2
    low_cutoff = 7; % 3 Hz
    high_cutoff = 35; % 35Hz
    b_fir = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    %figure
    %freqz(b_fir,1,512);
    
    % C3_filtered = filter(b_fir,1,C3);
    % C4_filtered = filter(b_fir,1,C4);
    % Cz_filtered = filter(b_fir,1,Cz);
    %
    % figure;
    % subplot(3,1,1)
    % plot(time,C3_filtered)
    % subplot(3,1,2)
    % plot(time,Cz_filtered)
    % subplot(3,1,3)
    % plot(time,C4_filtered)
    
    %% Spatial Filtering
    
    %%% 1ST : COMMON AVERAGE REFERENCE %%%%
    
    C3_prime = C3- (1/2)*(C4+Cz) ;
    C4_prime = C4- (1/2)*(C3+Cz) ;
    Cz_prime = Cz- (1/2)*(C4+C3) ;
    
    %%% 2ND : LAPLACIAN FILTER %%%
    d_C3 = [4 8]; %% Cz followed by C3, distance in cm
    w_hi = (1./d_C3)./(sum(1./d_C3));
    C3_lap = C3 - w_hi(1).*Cz - w_hi(2).*C4;
    
    %% Temporal Filtering
    
    C3_filtered = filter(b_fir_mu,1,C3_prime); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4_prime);
    Cz_filtered = filter(b_fir_mu,1,Cz_prime);
    
    figure
    subplot(3,1,1);
    plot(time,C3_filtered);
    title(' C3 After CAR & BDP');
    subplot(3,1,2);
    plot(time,Cz_filtered);
    title(' Cz After CAR & BDP');
    subplot(3,1,3);
    plot(time,C4_filtered);
    title(' C4 After CAR & BDP');
    
    %% CSP (Back to spatial filtering)
    
    %%% 3RD : Common spetcral patteren %%%. Notice that this is an adaptive
    %%% technique
    % Notice that we assume that X = [C3 Cz C4] (from left to right)
    
    R_ave_right = zeros(3,3);
    for i = 1:1:10
        load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Right_Trial_' num2str(i) '.mat']);
        R_normalized_current = spatial_cov_computation(Fs,C3,C4,Cz,time);
        R_ave_right = R_ave_right + R_normalized_current;
    end
    R_ave_right = R_ave_right./10;
    
    R_ave_left = zeros(3,3);
    for i = 1:1:10
        load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Left_Trial_' num2str(i) '.mat']);
        R_normalized_current = spatial_cov_computation(Fs,C3,C4,Cz,time);
        R_ave_left = R_ave_left + R_normalized_current;
    end
    R_ave_left = R_ave_left./10;
    R_tot = R_ave_right + R_ave_left;
    [U,D] = eig(R_tot);
    P = (inv(D))^(1/2)*(U.');
    %% Feature extraction using BandPower (or variance and a windows of 3s )
    %First Windowing around the moment
    t_move = 10;
    index_time = find(time ==t_move);
    wlen = 3*Fs;
    win = (hamming(wlen, 'periodic')).';
    
    C3_win = C3_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    C4_win = C4_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    Cz_win = Cz_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    var_tot = var(C3_win) + var(Cz_win) + var(C4_win); %% Used for normalisation ! 
    feature_vec = [var(C3_win); var(Cz_win); var(C4_win); 2*var_tot]./var_tot;
    yfit = trainedModel2.predictFcn(feature_vec(1:3)')
    feature_mat(i,:) = feature_vec.';
end


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
