function [R_normalized] = spatial_cov_computation(Fs,C3,C4,Cz,time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes her
t_move = 10;
index_time = find(time ==t_move);
wlen = 3*Fs;
win = (hamming(wlen, 'periodic')).';
C3_win = C3(index_time-Fs+1:index_time+2*Fs).*win;
C4_win = C4(index_time-Fs+1:index_time+2*Fs).*win;
Cz_win = Cz(index_time-Fs+1:index_time+2*Fs).*win;
X = [C3_win; Cz_win; C4_win];
R = X*(X.');
R_normalized = R./(trace(R));
end

