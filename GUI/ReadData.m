<<<<<<< HEAD
function [] = CreateDataFiles(file,number_class,Fs,n_channel)
=======
function [] = ReadData(file,number_class,Fs,n_channel,supindex)
>>>>>>> origin/Matthieu
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
Trecord = 3;
numFiles = round(size(data,1)/(Fs*Trecord));
secKept = 2*Fs;
t_vec = linspace(0,secKept/Fs,secKept);
<<<<<<< HEAD
for k = 1:1:numFiles
    fid  = fopen(['Newdata2/26_03_18_right' num2str(k) '.txt'],'w');
    for line = 1:1:secKept
        y = data((9 + (k-1)*(Fs*Trecord) + line),:);
        y = [y t_vec(line)];
        fprintf(fid,'%f\t %f\t %f\t %f\n',y);
    end 
    fclose(fid);
end

=======
% for k = 1:1:numFiles
%     fid  = fopen(['Newdata2/26_03_18_right' num2str(k) '.txt'],'w');
%     for line = 1:1:secKept
%         y = data((9 + (k-1)*(Fs*Trecord) + line),:);
%         y = [y t_vec(line)];
%         fprintf(fid,'%f\t %f\t %f\t %f\n',y);
%     end
%     fclose(fid);
% end
% number_data = size(data,1);
% t_vec = linspace(0,number_data/Fs,number_data);
% data = data.';
% C3 = data(1,:);
% Cz = data(2,:);
% C4 = data(3,:);
% time = t_vec;
% Rawdata = [C3; C4; Cz; time]
% save('Rawdata.mat','Rawdata');
% figure;
% subplot(3,1,1)
% plot(learning_set(:,1));
% subplot(3,1,2)
% plot(learning_set(:,2));
% subplot(3,1,3)
% plot(learning_set(:,3));
data = data.';
indexStr = 0;
for k = 1:1:numFiles
    indexStr = k + supindex;
    indx = 50 + (k-1)*(Fs*Trecord);
    C3 = data(1,indx+1:indx+secKept);
    Cz = data(2,indx+1:indx+secKept);
    C4 = data(3,indx+1:indx+secKept);
    %FC1 = data(4,indx+1:indx+secKept);
    %FC2 = data(6,indx+1:indx+secKept);
    CP2 = data(4,indx+1:indx+secKept);
    CP1 = data(5,indx+1:indx+secKept);
    time = t_vec;
    Rawdata = [C3; C4; Cz; time];
    save(['2504Data/25_04_18_Rest' num2str(indexStr) '.mat'],'C3','Cz','C4','CP2', 'CP1', 'time');
end
>>>>>>> origin/Matthieu
end

