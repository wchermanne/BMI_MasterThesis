clear all 
close all
clc

%% PARAMETERS
n_channel = 3;

n_hex = 7;
Gain = 24;
number_class = 2;
number_data = 1000;

% Initialization final matrix
%learning_set = zeros(number_class*number_data, n_channel+1);
learning_set = [];

% Left (1)
% DATA EXTRACTION

for i=4:4
    
file = fopen(['RealData/test' num2str(i) '.txt'],'r');

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

    Length = length(learning_set);
    learning_set(1+Length:Length+L,1:n_channel) = data;
    learning_set(1+Length:Length+L,n_channel+1) = 1;

end
figure;
subplot(3,1,1)
plot(learning_set(:,1));
subplot(3,1,2)
plot(learning_set(:,2));
subplot(3,1,3)
plot(learning_set(:,3));

%% Right (-1)
% DATA EXTRACTION 
for i=1:6
    
file = fopen(['29-03-17/right30s' num2str(i+1) '.txt'],'r');

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

    Length = length(learning_set);
    learning_set(1+Length:Length+L,1:3) = data;
    learning_set(1+Length:Length+L,4) = -1;

end

%% Rest (0)
% DATA EXTRACTION 
for i=1:6
    
file = fopen(['29-03-17/rest30s' num2str(i+1) '.txt'],'r');

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

    Length = length(learning_set);
    learning_set(1+Length:Length+L,1:3) = data;
    learning_set(1+Length:Length+L,4) = 0;

end

%% Save data set
save('learningset', 'learning_set');


