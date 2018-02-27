%%% This script is a complete BCI made of all the required stages :
% 1) FIR temporal filtering in the mu band (7-13Hz)
% 2) CSP Spatial Filter reducing the number of channel from three (C3,Cz,C4) to
% two (S1,S2)
% 3) Windowing & Framing stage
% 4) Feature Extraction based on log-normalized variance + Wavelet
% transform !!
% 5) Classification algorithm based on Linear SVM/ Linear LDA/ Fine KNN

%% TRAINING PART %%%%%%% Now 3 class are required ! Rigth, Left & Rest
side = 2;
if(side==1)
    text_side = 'Right';
else
    text_side = 'Left';
end
level = 5;
feature_mat = ones(20,3);
wavelet_feature_mat = ones(20,2*level+1);
wavelet_feature_mat_rest = ones(20,2*level+1);
feature_mat_rest = ones(20,3);
for i = 1:1:20
    %% Loading data
    load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/Data_for_CSP/EEG_Signals_' text_side '_Trial_' num2str(i) '.mat']);
    nyq_freq = Fs/2;
    %% Frequency Filtering & suuband creation for
    %%% One might use bandpass filtering from 3Hz (larger than EOG and eye blink artifacts)
    %%% up to 30-40 Hz (lower than power line interference @50 Hz)
    filter_order = 96;
    
    %%% 1ST : FIR Filter for mu rhythm --> pay attention to the delay = N/2
    low_cutoff = 7; % 3 Hz
    high_cutoff = 14; % 35Hz 14Hz
    b_fir_mu = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %%% 2ND : FIR Filter for beta rhythm --> pay attention to the delay = N/2
    low_cutoff = 15; % 3 Hz
    high_cutoff = 26; % 35Hz
    b_fir_beta = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
    
    %% Temporal Filtering
    
    C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4);
    Cz_filtered = filter(b_fir_mu,1,Cz);
    
    %% CSP (Back to spatial filtering)
    W_Csp = CSP_training(b_fir_mu);
    S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered];
    %% Feature extraction using BandPower (or variance and a windows of 3s )
    %First Windowing around the moment
    t_move = 10;
    index_time = find(time ==t_move);
    wlen = 3*Fs;
    win = (hamming(wlen, 'periodic')).';
    
    C3_win = C3_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    C4_win = C4_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    Cz_win = Cz_filtered(index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    %var_tot = var(C3_win) + var(Cz_win) + var(C4_win); %% Used for normalisation !
    %feature_vec = [var(C3_win); var(Cz_win); var(C4_win); 2*var_tot]./var_tot;
    %fit = trainedModel2.predictFcn(feature_vec(1:3)')
    %feature_mat(i,:) = feature_vec.';
    
    %%% Feature vector for Motor Task
    S1_win = S(1,index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    S2_win = S(2,index_time-Fs+1+filter_order/2:index_time+2*Fs+filter_order/2).*win;
    var_tot = var(S1_win) + var(S2_win); %% Used for normalisation !
    feature_vec = log([var(S1_win); var(S2_win)]./var_tot);
    wname = 'db1'; %% 'sym6' or 'db1'
    figure; subplot(2,1,1); plot(S1_win); subplot(2,1,2); plot(S2_win)
    
    [C_1,Level_1] = wavedec(S1_win,level,wname);
    [C_2,Level_2] = wavedec(S2_win,level,wname);
    wavelet_feature_vec = [];
    wavelet_var_tot = 0;
    for kl = 2:1:level+1
        current_length = Level_1(kl);
        Cplot =C_1(1:current_length);
        C_1 = C_1(current_length+1:end);
        current_length_2 = Level_2(kl);
        Cplot_2 =C_2(1:current_length_2);
        C_2 = C_2(current_length_2+1:end);
        wavelet_var_tot = wavelet_var_tot + var(Cplot) + var(Cplot_2);
        wavelet_feature_vec = [wavelet_feature_vec; (var(Cplot)) ; var(Cplot_2)];
    end
    
    wavelet_feature_mat(i,:) = [log(wavelet_feature_vec./wavelet_var_tot);side].';
    feature_mat(i,:) = [feature_vec;side].';
    
    S1_win = S(1,index_time-5*Fs+1+filter_order/2:index_time-2*Fs+filter_order/2).*win;
    S2_win = S(2,index_time-5*Fs+1+filter_order/2:index_time-2*Fs+filter_order/2).*win;
    var_tot = var(S1_win) + var(S2_win); %% Used for normalisation !
    feature_vec_rest = log([var(S1_win); var(S2_win)]./var_tot);
    
    figure; subplot(2,1,1); plot(S1_win); subplot(2,1,2); plot(S2_win)
    [C_1,Level_1] = wavedec(S2_win,level,wname);
    [C_2,Level_2] = wavedec(S2_win,level,wname);
    wavelet_feature_vec = [];
    wavelet_var_tot = 0;
    
    figure;
    for kl = 2:1:level+1
        current_length = Level_1(kl);
        Cplot =C_1(1:current_length);
        C_1 = C_1(current_length+1:end);
        current_length_2 = Level_2(kl);
        Cplot_2 =C_2(1:current_length_2);
        C_2 = C_2(current_length_2+1:end);
        wavelet_var_tot = wavelet_var_tot + var(Cplot) + var(Cplot_2);
        wavelet_feature_vec = [wavelet_feature_vec; (var(Cplot)) ; var(Cplot_2)];
    end
    wavelet_feature_mat_rest(i,:) = [log(wavelet_feature_vec./wavelet_var_tot);3].';
    feature_mat_rest(i,:) = [feature_vec_rest;3].';
    
end
close all

%%% First trial of machine learning prediction






%% TESTING PART %%%%%%%%%%%%%
load('KNN.mat')
load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/Data_for_CSP/EEG_Signals_Left_Trial_6.mat']);

%%% Temporal Filtering

C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
C4_filtered = filter(b_fir_mu,1,C4);
Cz_filtered = filter(b_fir_mu,1,Cz);

%%% CSP
S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered];
%%% Wavelet Transform


wname = 'sym1'; %% 'sym6' or 'db1'
level = 6;

[C_1,Level_1] = wavedec(S(1,:),level,wname);
figure;
for kl = 1:1:level+1
    subplot(level+1,1,kl)
    current_length = Level_1(kl);
    Cplot =C_1(1:current_length)
    C_1 = C_1(current_length+1:end);
    plot(Cplot)
    title(['cD' num2str(kl)]);
end


[C_2,Level_2] = wavedec(S(2,:),level,wname);
figure;
for kl = 1:1:level+1
    subplot(level+1,1,kl)
    current_length = Level_2(kl);
    Cplot =C_2(1:current_length)
    C_2 = C_2(current_length+1:end);
    plot(Cplot)
    title(['cD' num2str(kl)]);
end

figure; subplot(2,1,1); plot(S(1,:)); subplot(2,1,2); plot(S(2,:));
%%% Windowing & Framing
WindowedSig = windowing_framing(Fs,time,S,3,1000);
S1 =  WindowedSig.Channel(1).X;
S2 = WindowedSig.Channel(2).X;

for k = 1:1:size(WindowedSig.Channel(1).X,2)
    S1_win = S1(:,k);
    S2_win = S2(:,k);
    var_tot = var(S1_win) + var(S2_win); %% Used for normalisation !
    feature_vec = log([var(S1_win); var(S2_win)]./var_tot);
    fit = trainedModelKNN.predictFcn(feature_vec');
    if(fit == 1)
        char = 'Left Movement Detected!'
    elseif(fit==2)
        char = 'Right Movement Detected!'
    end
end
