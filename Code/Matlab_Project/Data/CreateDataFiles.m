function [] = CreateDataFiles(file,number_class,Fs,n_channel)
%% PARAMETERS
n_hex = 7;
Gain = 24;

% Initialization final matrix
learning_set = [];

% Left (1)
% DATA EXTRACTION


file = fopen([file] ,'r');

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
Trecord = 4;
numFiles = round(size(data,1)/(Fs*Trecord));
secKept = 3.5*Fs;
t_vec = linspace(0,secKept/Fs,secKept);

data = data.';
mydate=date;
for k = 1:1:numFiles
    indx = 9 + (k-1)*(Fs*Trecord);
    C3 = data(1,indx+1:indx+secKept);
    Cz = data(2,indx+1:indx+secKept);
    C4 = data(3,indx+1:indx+secKept);
    time = t_vec;
    Rawdata = [C3; C4; Cz; time];
    save([mydate '_right' num2str(k) '.mat'],'C3','Cz','C4', 'time');
end
end

