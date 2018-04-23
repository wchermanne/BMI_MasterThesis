%%% This script is a complete BCI made of all the required stages :
% 1) FIR temporal filtering in the mu band (7-13Hz)
% 2) CSP Spatial Filter reducing the number of channel from three (C3,Cz,C4) to
% two (S1,S2)
% 3) Windowing & Framing stage
% 4) Feature Extraction based on log-normalized variance + Wavelet
% transform !!
% 5) Classification algorithm based on Linear SVM/ Linear LDA/ Fine KNN

%% TRAINING PART %%%%%%% Now 3 class are required ! Rigth, Left & Rest
clear all
dwtmode('sym');
TrainingSet =[];
TrainingSet2 = [];
predictionMat = [];
for side = 1:1:3 
    if(side==1)
        text_side = 'Right';
    elseif(side ==2)
        text_side = 'Left';
    else
        text_side = 'Rest';
    end
    Fs = 250;
    level = 4;
    feature_mat = ones(20,(2*1)+1);
    wavelet_feature_mat = ones(20,(6*(level))*1+1);
    Fs = Fs/2;
    nyq_freq = Fs/2;
    for i = 1:1:40
        %% Loading data
        feature_mat_bands = [];
        wavelet_feature_vec_bands = [];
        load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/RealData7channels/13_04_18_' text_side  num2str(i) '.mat']);
        C3 = downsample(C3,2);
        C4 = downsample(C4,2);
        Cz = downsample(Cz,2);
        CP1 = downsample(CP1,2);
        CP2 = downsample(CP2,2);
        FC1 = downsample(FC1,2);
        FC2 = downsample(FC2,2);
        time = downsample(time,2);
        for freq =1:1:1
            %% Frequency Filtering & suuband creation for
            if (freq ==1)
                low_cutoff = 7;%7; % 3 Hz
                high_cutoff = 14;%14; % 35Hz 14Hz
            else
                low_cutoff = 18; % 3 Hz
                high_cutoff = 26; % 35Hz 14Hz
            end
            
            %%% One might use bandpass filtering from 3Hz (larger than EOG and eye blink artifacts)
            %%% up to 30-40 Hz (lower than power line interference @50 Hz)
            filter_order = 20;
            
            %%% 1ST : FIR Filter for mu rhythm --> pay attention to the delay = N/2
            b_fir_mu = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);
            [b_iir_mu,a_iir_mu] = butter(16,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
            
            %%% 2ND : FIR Filter for beta rhythm --> pay attention to the delay = N/2
            b_fir_beta = fir1(filter_order,[low_cutoff./nyq_freq,high_cutoff./nyq_freq]);% Hamming Window-based FIR filter design (from 0->1 with 1 correpons to nyquist frequency)
            %% Notch Filtering
            wo = 50/(Fs/2);  bw = wo/35;
            [b,a] = iirnotch(wo,bw);
%             C3_filtered = filter(b,a,C3); % C3_prime but let's modify it
%             C4_filtered = filter(b,a,C4);
%             Cz_filtered = filter(b,a,Cz);
            C3_filtered = C3; % C3_prime but let's modify it
            C4_filtered = C4;
            Cz_filtered = Cz;
            CP1_filtered = CP1;
            CP2_filtered = CP2;
            FC1_filtered = FC1;
            FC2_filtered = FC2;
            %% Temporal Filtering
            
%             C3_filtered = filter(b_fir_mu,1,C3_filtered); % C3_prime but let's modify it
%             C4_filtered = filter(b_fir_mu,1,C4_filtered);
%             Cz_filtered = filter(b_fir_mu,1,Cz_filtered);
            C3_filtered = filter(b_fir_mu,1,C3_filtered); % C3_prime but let's modify it
            C4_filtered = filter(b_fir_mu,1,C4_filtered);
            Cz_filtered = filter(b_fir_mu,1,Cz_filtered);
            CP1_filtered = filter(b_fir_mu,1,CP1_filtered);
            CP2_filtered = filter(b_fir_mu,1,CP2_filtered);
            FC1_filtered = filter(b_fir_mu,1,FC1_filtered);
            FC2_filtered = filter(b_fir_mu,1,FC2_filtered);
            
            %% CSP (Back to spatial filtering)
            W_Csp = CSP_training_Real_Data_Downsampling_7_channles(b_fir_mu);
            S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered; FC1_filtered; FC2_filtered; CP1_filtered; CP2_filtered]; % [C3_win; Cz_win; C4_win; FC1; FC2; CP1; CP2];
            %S = fastica([C3_filtered; Cz_filtered; C4_filtered],'NumOfIC',2);
            %% Feature extraction using BandPower (or variance and a windows of 3s )
            %First Windowing around the moment
            t_move = 10;
            index_time = find(time==t_move);
            wlen = 3*Fs;
            win = (hamming(wlen, 'periodic')).';
            
            %%% Feature vector for Motor Task
            S1_win = S(1,:);
            S2_win = S(2,:);
            var_tot = var(S1_win) + var(S2_win); %% Used for normalisation !
            feature_vec = log([var(S1_win); var(S2_win)]./var_tot);
            wname = 'db1'; %% 'sym6' or 'db3'
            %figure; subplot(2,1,1); plot(S1_win); subplot(2,1,2); plot(S2_win)
            %%%%%%%% /! \ %%%%%%% CHANGE MODE IF NECESSARY
            
            %dwtmode('sym');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            %LO_D = [0.0352262918821007 -0.0854412738822415 -0.135011020010391 0.459877502119331 0.806891509313339 0.332670552950957];
            %HI_D = [-0.332670552950957 0.806891509313339 -0.459877502119331 -0.135011020010391 0.0854412738822415 0.0352262918821007];
            [C_1,Level_1] = wavedec(S1_win,level,wname);
            [C_2,Level_2] = wavedec(S2_win,level,wname);
            wavelet_feature_vec = [];
            wavelet_feature_vec_mean = [];
            wavelet_feature_vec_min = [];
            wavelet_feature_vec_max = [];
            wavelet_var_tot = 0;
            wavelet_mean_tot = 0;
            wavelet_min_tot = 0;
            wavelet_max_tot = 0;
            C1_copy = C_1 ;
            C2_copy = C_2;
            C_1 = C_1(Level_1(1)+1:end);
            C_2 = C_2(Level_2(1)+1:end);
            for kl = 2:1:level+1 %% Keep only decomposition coeffeicients 
                current_length = Level_1(kl);
                Cplot =C_1(1:current_length);
                C_1 = C_1(current_length+1:end);
                current_length_2 = Level_2(kl);
                Cplot_2 =C_2(1:current_length_2);
                C_2 = C_2(current_length_2+1:end);
                wavelet_var_tot = wavelet_var_tot + var(Cplot) + var(Cplot_2);
                wavelet_mean_tot = wavelet_mean_tot + mean(Cplot) + mean(Cplot_2);
                wavelet_min_tot =   wavelet_min_tot + min(Cplot) + min(Cplot_2);
                wavelet_max_tot =   wavelet_max_tot + max(Cplot) + max(Cplot_2);
                wavelet_feature_vec = [wavelet_feature_vec; (var(Cplot)) ; var(Cplot_2)];
                wavelet_feature_vec_mean = [wavelet_feature_vec_mean; mean(Cplot); mean(Cplot_2)];
                wavelet_feature_vec_min = [wavelet_feature_vec_min; min(Cplot); min(Cplot_2)];
                wavelet_feature_vec_max = [wavelet_feature_vec_max; max(Cplot); max(Cplot_2)];
            end
            
            %wavelet_feature_mat(i,:) = [log(wavelet_feature_vec./wavelet_var_tot);side].';
            %feature_mat(i,:) = [feature_vec;side].';
            feature_mat_bands = [feature_mat_bands; feature_vec];
            wavelet_feature_vec_bands = [wavelet_feature_vec_bands ; log(wavelet_feature_vec./wavelet_var_tot); (wavelet_feature_vec_max./wavelet_max_tot); (wavelet_feature_vec_min./wavelet_min_tot);  ];
        end
        feature_mat(i,:) = [feature_mat_bands;side].';
        wavelet_feature_mat(i,:) = [wavelet_feature_vec_bands;side].';
    end
    TrainingSet = [feature_mat; TrainingSet];
    TrainingSet2 = [wavelet_feature_mat; TrainingSet2];
    %[y1] = myNeuralNetworkFunctionReduced(wavelet_feature_mat(:,1:end-1));
    %predictionMat = [y1; predictionMat];
end
TrainngSetNN = TrainingSet2(:,1:end-1);
TrainingSetTarget = zeros(120,3);
TrainingSetTarget(1:40,3) = ones(40,1);
TrainingSetTarget(41:80,2) = ones(40,1);
TrainingSetTarget(81:120,1) = ones(40,1);
close all

%%% First trial of machine learning prediction






%% TESTING PART %%%%%%%%%%%%%
load('KNN.mat')
load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/Data_for_CSP/EEG_Signals_Right_Trial_18.mat']);

%%% Temporal Filtering

C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
C4_filtered = filter(b_fir_mu,1,C4);
Cz_filtered = filter(b_fir_mu,1,Cz);

%%% CSP
S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered];
%%% Wavelet Transform
figure;
subplot(2,1,1)
plot(S(1,:));
subplot(2,1,2)
plot(S(2,:));
% wname = 'sym1'; %% 'sym6' or 'db1'
% level = 6;
%
% [C_1,Level_1] = wavedec(S(1,:),level,wname);
% figure;
% for kl = 1:1:level+1
%     subplot(level+1,1,kl)
%     current_length = Level_1(kl);
%     Cplot =C_1(1:current_length)
%     C_1 = C_1(current_length+1:end);
%     plot(Cplot)
%     title(['cD' num2str(kl)]);
% end
%
%
% [C_2,Level_2] = wavedec(S(2,:),level,wname);
% figure;
% for kl = 1:1:level+1
%     subplot(level+1,1,kl)
%     current_length = Level_2(kl);
%     Cplot =C_2(1:current_length)
%     C_2 = C_2(current_length+1:end);
%     plot(Cplot)
%     title(['cD' num2str(kl)]);
% end
%
% figure; subplot(2,1,1); plot(S(1,:)); subplot(2,1,2); plot(S(2,:));

%%% Windowing & Framing
WindowedSig = windowing_framing(Fs,time,S,3,1000);
S1 =  WindowedSig.Channel(1).X;
S2 = WindowedSig.Channel(2).X;
timeS = WindowedSig.time;
for k = 1:1:size(WindowedSig.Channel(1).X,2)
    S1_win = S1(:,k);
    S2_win = S2(:,k);
    timeSwin = timeS(:,k);
    wname = 'db1';
    feature_vec = FvWavelets([S1_win S2_win],timeSwin,Fs,level,wname);
    %     [C_1,Level_1] = wavedec(S1_win,level,wname);
    %     [C_2,Level_2] = wavedec(S2_win,level,wname);
    %     wavelet_feature_vec = [];
    %     wavelet_var_tot = 0;
    %     for kl = 2:1:level+1
    %         current_length = Level_1(kl);
    %         Cplot =C_1(1:current_length);
    %         C_1 = C_1(current_length+1:end);
    %         current_length_2 = Level_2(kl);
    %         Cplot_2 =C_2(1:current_length_2);
    %         C_2 = C_2(current_length_2+1:end);
    %         wavelet_var_tot = wavelet_var_tot + var(Cplot) + var(Cplot_2);
    %         wavelet_feature_vec = [wavelet_feature_vec; (var(Cplot)) ; var(Cplot_2)];
    %     end
    %     feature_vec = log(wavelet_feature_vec./wavelet_var_tot);
    
    
    fit = trainedModel.predictFcn(feature_vec.');
    fit_2 = FtNaiveKNN(trainedModel,feature_vec);
    if (fit_2 ~= fit)
        char = 'error';
    end
    
    if(fit_2 == 1)
        char = 'Right Movement Detected!'
    elseif(fit_2==2)
        char = 'Left Movement Detected!'
    elseif(fit_2 ==3)
        char = 'Rest'
    end
end

