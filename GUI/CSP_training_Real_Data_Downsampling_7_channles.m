function [W_Csp] = CSP_training_Real_Data_Downsampling_7_channles(b_fir_mu)
%This function performs the training of the CSP filter W
%%% Common spetcral patteren %%%. Notice that this is an adaptive
%%% technique
% Notice that we assume that X = [C3 Cz C4] (from left to right)

R_ave_right = zeros(7,7);
for k = 1:1:40
    load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/RealData7channels/13_04_18_Right'  num2str(k) '.mat']);
    C3 = downsample(C3,2);
    C4 = downsample(C4,2);
    Cz = downsample(Cz,2);
    CP1 = downsample(CP1,2);
    CP2 = downsample(CP2,2);
    FC1 = downsample(FC1,2);
    FC2 = downsample(FC2,2);
    C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4);
    Cz_filtered = filter(b_fir_mu,1,Cz);
    FC1_filtered = filter(b_fir_mu,1,FC1);
    FC2_filtered = filter(b_fir_mu,1,FC2);
    CP1_filtered = filter(b_fir_mu,1,CP1);
    CP2_filtered = filter(b_fir_mu,1,CP2);
    R_normalized_current = spatial_cov_computation_Real_Data_7_channels(C3_filtered,C4_filtered,Cz_filtered,  FC1_filtered ,  FC2_filtered,   CP1_filtered,   CP2_filtered);
    R_ave_right = R_ave_right + R_normalized_current;
end
R_ave_right = R_ave_right./40;

R_ave_left = zeros(7,7);
for k = 1:1:40
    load(['/Users/matthieu/GitHub/BMI_MasterThesis/GUI/RealData7channels/13_04_18_Left'  num2str(k) '.mat']);
    C3 = downsample(C3,2);
    C4 = downsample(C4,2);
    Cz = downsample(Cz,2);
    CP1 = downsample(CP1,2);
    CP2 = downsample(CP2,2);
    FC1 = downsample(FC1,2);
    FC2 = downsample(FC2,2);
    C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4);
    Cz_filtered = filter(b_fir_mu,1,Cz);
    FC1_filtered = filter(b_fir_mu,1,FC1);
    FC2_filtered = filter(b_fir_mu,1,FC2);
    CP1_filtered = filter(b_fir_mu,1,CP1);
    CP2_filtered = filter(b_fir_mu,1,CP2);
    R_normalized_current = spatial_cov_computation_Real_Data_7_channels(C3_filtered,C4_filtered,Cz_filtered,  FC1_filtered ,  FC2_filtered,   CP1_filtered,   CP2_filtered);
    R_ave_left = R_ave_left + R_normalized_current;
end
R_ave_left = R_ave_left./40;
R_tot = R_ave_right + R_ave_left;
[U,D] = eig(R_tot);
P = (inv(D))^(1/2)*(U.');
Sigma_1_hat = P*R_ave_right*(P.');
Sigma_2_hat = P*R_ave_left*(P.'); % The sum of the S1 & S2 = I
[V,Gamma] = eig(Sigma_1_hat);
W = (P.')*V;
%%%% Then we can select first and last vector. They have max variance
%%%% for class 1 & 2 respectively
W_Csp = [W(:,1) W(:,7)];
%S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered];
end





