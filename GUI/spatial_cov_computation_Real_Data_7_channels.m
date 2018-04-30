function [R_normalized] = spatial_cov_computation_Real_Data_7_channels(C3,C4,Cz,FC1,FC2,CP1,CP2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes her
% t_move = 10;
% index_time = find(time ==t_move);
% wlen = 3*Fs;
% win = (hamming(wlen, 'periodic')).';
C3_win = C3;
C4_win = C4;
Cz_win = Cz;
%X = [C3_win; Cz_win; C4_win; FC1; FC2; CP1; CP2];
X = [C3_win; Cz_win; C4_win; CP2; CP1];
R = X*(X.');
R_normalized = R./(trace(R));
end

