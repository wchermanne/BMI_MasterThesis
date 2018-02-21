function [C3,C4,Cz,time] = eeg_generator(Fs,t_end,SNR,movement, time_movement,eye_blink,power_line_check)
%%% Signal generation for test (continuous sensorimotor rythm ) using 1st Order MPA
%%% LEFT = 2
%%% RIGHT =1
freq = linspace(0,41,42);
amp = [0.5 4/8 7/8 4/8 2.5 0.3 0.15 0.3 3 5 6 8 5 4 2 1.5 1.2 1 0.9 0.6 0.4 0.2 0.1 0.1 0.05 0.03 0.02 0.01 0.008 0.004 0.002 0.001 0.0008 0.0005 0.0003 0.0001 0.0001 0.00008 0.00008 0.000053 0.000033 0.00002]*10;
if(power_line_check ==1)
    freq = [freq 50]';
    amp = [amp 200];
end
amp_4 = amp;
amp_z = amp;

phi = linspace(pi/4,pi/2,length(amp)) + pi./10*rand(1,length(amp)) - pi./20;
t = 0:1/Fs:t_end;
sin_mat_C3 = zeros(1,length(t));
sin_mat_C4 = zeros(1,length(t));
sin_mat_Cz = zeros(1,length(t));
logbook_gamma = zeros(1,length(t));
logbook_amp =  zeros(1,length(t));

%%%% Lower mu rhythm %%%%%
mu_ERD_duration = 0.5 + 0.1*rand;

%%%% Lower mu rhythm %%%%%
upper_mu_ERS = (time_movement./2)*(0.6 + 0.2*rand); %% between 60% and 100%
upper_mu_ERD = (time_movement./2)*(0.8 + 0.4*rand); %% Between 80% and 100 %

%%%% Characteristics of the movement %%%%
movement_duration = 1;

%%% Taking eye blink effect into account
eye_blink_occurence = eye_blink;
time_eye_blink = ones(1,eye_blink_occurence);
eye_blink_rising_duration = 0.08 + 0.01*rand;
eye_blink_steady_duration = ones(1,eye_blink_occurence);
for k =1:1:eye_blink_occurence
    time_eye_blink(k) = 1+(t_end./2)*k*rand;
    eye_blink_steady_duration(k) =0.2 + 0.2*rand;
end

%%% Gamma ERD prior to a movement
gamma_ERD_duration = 0.3 + 0.1*rand;

for l = 1:1:length(t)
    %%% Taking Event-related desynchronization and resync for lower mu
    %%% rhythm. That happens for any motor task
    if (t(l) <= time_movement-mu_ERD_duration && t(l) >= time_movement-2*mu_ERD_duration)
        amp(8:10) = 0.98.*amp(8:10);
        amp_4(8:10) = 0.98.*amp_4(8:10);
        amp_z(8:10) = 0.98.*amp_z(8:10);
    end
    if (t(l) <= time_movement && t(l) >= time_movement-mu_ERD_duration)
        amp(8:10) = 1.*amp(8:10)./0.98;
        amp_4(8:10) = 1.*amp_4(8:10)./0.98;
        amp_z(8:10) = 1.*amp_z(8:10)./0.98;
    end
    
    %%% Taking ERS for upper mu (ipsilateral side increases)
    %%% rhythm.
    ERS_mu_factor = 1.002;
    ERS_mu_factor_z = 1.0005;
    if (t(l) >= time_movement && t(l) <= time_movement + upper_mu_ERS/2)
        if(movement == 2)
            amp(11:13) = ERS_mu_factor .*amp(11:13);
            amp_z(11:13) = ERS_mu_factor_z.*amp_z(11:13);
        elseif(movement == 1)
            amp_4(11:13) = ERS_mu_factor .*amp_4(11:13);
            amp_z(11:13) = ERS_mu_factor_z.*amp_z(11:13);
        end
    end
    if (t(l) >= time_movement + upper_mu_ERS/2 && t(l) <= time_movement + upper_mu_ERS)
        if(movement == 2)
            amp(11:13) = amp(11:13)./(ERS_mu_factor );
            amp_z(11:13) = amp_z(11:13)./(ERS_mu_factor_z);
        elseif(movement == 1)
            amp_4(11:13) = amp_4(11:13)./(ERS_mu_factor );
            amp_z(11:13) = amp_z(11:13)./(ERS_mu_factor_z);
        end
    end
    
        %%% Taking ERD for upper mu
    %%% rhythm.
    if (t(l) >= time_movement && t(l) <= time_movement + upper_mu_ERD/2)
        if(movement == 2)
            amp_4(11:13) = 0.985.*amp_4(11:13);
        elseif(movement == 1)
            amp(11:13) = 0.985.*amp(11:13);
        end
    end
    if (t(l) >= time_movement + upper_mu_ERD/2 && t(l) <= time_movement + upper_mu_ERD)
        if(movement == 2)
            amp_4(11:13) = amp_4(11:13)./(0.985);
        elseif(movement == 1)
            amp(11:13) = amp(11:13)./(0.985);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Taking Beta rhythm into account
    if (t(l) >= time_movement && t(l) <= time_movement+movement_duration)
        if(movement == 2)
            amp(18:26) = 1.005.*amp(18:26);
            amp_4(18:26) = 1.004.*amp_4(18:26);
            amp_z(18:26) = 1.0045.*amp_z(18:26);
        elseif (movement == 1)
            amp(18:26) = 1.004.*amp(18:26);
            amp_4(18:26) = 1.005.*amp_4(18:26);
            amp_z(18:26) = 1.0045.*amp_z(18:26);
        end
    end
    if (t(l) >= time_movement+movement_duration && t(l) <= time_movement+2*movement_duration)
        if(movement == 2)
            amp(18:26) = amp(18:26)./(1.005);
            amp_4(18:26) = amp_4(18:26)./(1.004);
            amp_z(18:26) = amp_z(18:26)./(1.0045);
        elseif (movement == 1)
            amp(18:26) = amp(18:26)./(1.004);
            amp_4(18:26) = amp_4(18:26)./(1.005);
            amp_z(18:26) = amp_z(18:26)./(1.0045);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Taking Gamma rhythm into account
    if (t(l) <= time_movement && t(l) >= time_movement-1*gamma_ERD_duration)
        amp(39:41) = 1.11.*amp(39:41);
        amp_4(39:41) = 1.11.*amp_4(39:41);
        amp_z(39:41) = 1.11.*amp_z(39:41);
    end
    if (t(l) <= time_movement+gamma_ERD_duration && t(l) >= time_movement)
        amp(39:41) = amp(39:41)./(1.11);
        amp_4(39:41) = amp_4(39:41)./(1.11);
        amp_z(39:41) = amp_z(39:41)./(1.11);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Taking Eye Blink
    if(isempty(time_eye_blink) ==1)
        
    else
        for k =1:1:eye_blink_occurence
            if(t(l) >= time_eye_blink(k) && t(l) <=time_eye_blink(k)+eye_blink_rising_duration)
                amp(1:2) = 1.3.*amp(1:2);
                amp_4(1:2) = 1.3.*amp_4(1:2);
                amp_z(1:2) = 1.3.*amp_z(1:2);
            end
            if(t(l) >= time_eye_blink(k) + eye_blink_rising_duration +eye_blink_steady_duration && t(l) <=time_eye_blink(k)+2*eye_blink_rising_duration+eye_blink_steady_duration)
                amp(1:2) = 1.*amp(1:2)./(1.3);
                amp_4(1:2) = 1.*amp_4(1:2)./(1.3);
                amp_z(1:2) = 1.*amp_z(1:2)./(1.3);
            end
        end
        
    end
    logbook_amp(l) = amp(20);
    sin_mat_temp = 0;
    sin_mat_temp_4 = 0;
    sin_mat_temp_z = 0;
    for k = 1:1:length(freq)
        sin_mat_temp = sin_mat_temp + amp(k)*sin(2*pi*freq(k)*t(l)+phi(k));
        sin_mat_temp_4 = sin_mat_temp_4 + amp_4(k)*sin(2*pi*freq(k)*t(l)+phi(k));
        sin_mat_temp_z = sin_mat_temp_z + amp_z(k)*sin(2*pi*freq(k)*t(l)+phi(k));
    end
    sin_mat_C3(l) = sin_mat_temp;
    sin_mat_C4(l) = sin_mat_temp_4;
    sin_mat_Cz(l) = sin_mat_temp_z;
    gamma = 0.995 + 0.01*rand;
    logbook_gamma(l) = gamma;
    amp = gamma*amp; + 4*rand(1,length(amp)) - 2;
    amp_4 = gamma*amp_4; + 4*rand(1,length(amp_4)) - 2;
    amp_z = gamma*amp_z; + 4*rand(1,length(amp_z)) - 2;
    if(power_line_check ==1)
        amp(end) = 200;
        amp_4(end) = 200;
        amp_z(end) = 200;
    end
end
% figure; plot(logbook_amp)
%%% Add noise contribution
Y = awgn(sin_mat_C3,SNR,'measured');
C3= Y;
C4 = awgn(sin_mat_C4,SNR,'measured');
Cz = awgn(sin_mat_Cz,SNR,'measured');

%% Save the signals and parameters
% save('channels.mat','C3','C4','Cz')
% 
% parameters=[];
% parameters.duration=duration;
% parameters.fsample=Fs;
% parameters.timestart=time_movement;
% parameters.movement=movement;
% parameters.eye_blink_check=eye_blink;
% parameters.power_line_check=power_line_check;
% save('parameters.mat', 'parameters')
% 

%% Plots and figures
time = t;
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
wlen = 2*Fs;
hop = Fs;
[stft_tot, freq_vec, time_vec] = stft(Y, wlen, hop, dft_size, Fs);
figure;
surf(time_vec,freq_vec',mag2db(abs(stft_tot)));

figure;
spectrogram(C3,512,256,dft_size,Fs,'psd');
title('C3 STFT');

figure;
spectrogram(C4,512,256,dft_size,Fs,'psd');
title('C4 STFT')
% nyq_freq = Fs./2; %% Half the sampling rate; nyquist frequency
% wlen = 2*Fs;
% hop = Fs;
% [stft_tot, freq_vec, time_vec] = stft(Y, wlen, hop, dft_size, Fs);
% figure;
% surf(time_vec,freq_vec',mag2db(abs(stft_tot)));


figure;
spectrogram(C3,512,256,dft_size,Fs,'psd');
title('C3 STFT');

figure;
spectrogram(C4,512,256,dft_size,Fs,'psd');
title('C4 STFT')
end

