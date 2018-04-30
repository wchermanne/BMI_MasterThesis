function [R_normalized] = spatial_cov_computation_Real_Data(C3,C4,Cz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes her
% t_move = 10;
% index_time = find(time ==t_move);
% wlen = 3*Fs;
% win = (hamming(wlen, 'periodic')).';
C3_win = C3;
C4_win = C4;
Cz_win = Cz; 
X = [C3_win; Cz_win; C4_win];
R = X*(X.');
R_normalized = R./(trace(R));
end

