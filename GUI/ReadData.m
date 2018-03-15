function [C3,Cz,C4,time] = ReadData(file,number_class,Fs,n_channel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% PARAMETERS
n_hex = 7;
Gain = 24;

% Initialization final matrix
learning_set = [];

% Left (1)
% DATA EXTRACTION


file = fopen(['RealData/' file] ,'r');

nextline = '';
str='';
while ischar(nextline)
    nextline = fgetl(file);
    if ( ischar(nextline) && (length(nextline)==n_channel*(n_hex+1)) )
        str = [str;nextline];
    end
end
fclose(file);

% DATA CONVERSION
L = length(str);
Vref = 4500;
LSB = Vref/(2^23-1)/Gain;

data = zeros(L,n_channel);

tmp = sfi([],25,0);
tmp_char = '';
for i=1:L
    for j=0:n_channel-1
        tmp_char = str(i,(j*(n_hex+1)+1):(j+1)*n_hex+j);
        tmp.hex = tmp_char;
        data(i,j+1) = tmp.data;
    end
end
data = data.*LSB;

Length = length(learning_set);
learning_set(1+Length:Length+L,1:n_channel) = data;
learning_set(1+Length:Length+L,n_channel+1) = 1;

number_data = size(data,1);
t_vec = linspace(0,number_data/Fs,number_data);
data = data.';
C3 = data(1,:);
Cz = data(2,:);
C4 = data(3,:);
time = t_vec;
Rawdata = [C3; C4; Cz; time]
save('Rawdara.mat','Rawdata');
% figure;
% subplot(3,1,1)
% plot(learning_set(:,1));
% subplot(3,1,2)
% plot(learning_set(:,2));
% subplot(3,1,3)
% plot(learning_set(:,3));
end

