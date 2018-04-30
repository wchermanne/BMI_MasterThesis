function [W_Csp] = CSP_training(b_fir_mu)
%This function performs the training of the CSP filter W
%%% Common spetcral patteren %%%. Notice that this is an adaptive
%%% technique
% Notice that we assume that X = [C3 Cz C4] (from left to right)

R_ave_right = zeros(3,3);
for k = 1:1:10
    load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Right_Trial_' num2str(k) '.mat']);
    C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4);
    Cz_filtered = filter(b_fir_mu,1,Cz);
    R_normalized_current = spatial_cov_computation(Fs,C3_filtered,C4_filtered,Cz_filtered,time);
    R_ave_right = R_ave_right + R_normalized_current;
end
R_ave_right = R_ave_right./10;

R_ave_left = zeros(3,3);
for k = 1:1:10
    load(['/Users/matthieu/GitHub/MasterThesis_BCI/GUI/Data_for_CSP/EEG_Signals_Left_Trial_' num2str(k) '.mat']);
    C3_filtered = filter(b_fir_mu,1,C3); % C3_prime but let's modify it
    C4_filtered = filter(b_fir_mu,1,C4);
    Cz_filtered = filter(b_fir_mu,1,Cz);
    R_normalized_current = spatial_cov_computation(Fs,C3_filtered,C4_filtered,Cz_filtered,time);
    R_ave_left = R_ave_left + R_normalized_current;
end
R_ave_left = R_ave_left./10;
R_tot = R_ave_right + R_ave_left;
[U,D] = eig(R_tot);
P = (inv(D))^(1/2)*(U.');
Sigma_1_hat = P*R_ave_right*(P.');
Sigma_2_hat = P*R_ave_left*(P.'); % The sum of the S1 & S2 = I 
[V,Gamma] = eig(Sigma_1_hat);
W = (P.')*V;
%%%% Then we can select first and last vector. They have max variance
%%%% for class 1 & 2 respectively
W_Csp = [W(:,1) W(:,3)];
%S = (W_Csp.')*[C3_filtered; Cz_filtered; C4_filtered];
end

