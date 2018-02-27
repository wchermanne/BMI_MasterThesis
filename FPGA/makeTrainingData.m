clear all 
close all

%% PARAMETERS
n_channel = 3;

n_hex = 7;
Gain = 24;
number_class = 3;
number_data = 15000;

%% Initialization final matrix
learning_set = zeros(number_class*number_data, n_channel+1);

%% Left (1)
% DATA EXTRACTION
file = fopen('Acquisition_Left_60s.txt','r');

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
Fs = 250;
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

learning_set(1:number_data,1:3) = data;
learning_set(1:number_data,4) = 1;

%% Right (-1)
% DATA EXTRACTION 
file = fopen('Acquisition_Right_60s.txt','r');

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
Fs = 250;
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

learning_set(number_data+1:2*number_data,1:3) = data;
learning_set(number_data+1:2*number_data,4) = -1;

%% Rest (0)
% DATA EXTRACTION 
file = fopen('Acquisition_Rest_60s.txt','r');

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
Fs = 250;
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

learning_set(2*number_data+1:3*number_data,1:3) = data;
learning_set(2*number_data+1:3*number_data,4) = 0;

%% Save data set
save('learningset', 'learning_set');


