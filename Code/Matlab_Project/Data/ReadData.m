function [data,time] = ReadData(file,number_class,Fs,n_channel,figure)

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

number_data = size(data,1)-250; % Modified to limit to 10s
t_vec = linspace(0,number_data/Fs,number_data);
data = data.';
data=data(:,1:2500);
C3 = data(1,1:2500); % Modified % Modified to limit to 10s
Cz = data(2,1:2500); % Modified % Modified to limit to 10s
C4 = data(3,1:2500);
if(n_channel==11)
    Fc5 = data(4,1:2500); % Modified % Modified to limit to 10s
    Fc1 = data(5,1:2500); % Modified % Modified to limit to 10s
    Fc2 = data(6,1:2500);
    Fc6 = data(7,1:2500);
    Cp6 = data(8,1:2500);
    Cp2 = data(9,1:2500);
    Cp1 = data(10,1:2500);
    Cp5 = data(11,1:2500);
end
assignin('base','data',data);
assignin('base','time',t_vec);

time = t_vec;

if(figure==1 && n_channel==3)
figure;
subplot(3,1,1)
plot(t_vec,C3);
title('C3')
xlabel('time [s]')
subplot(3,1,2)
plot(t_vec,Cz);
title('Cz')
xlabel('time [s]')
subplot(3,1,3)
plot(t_vec,C4);
title('C4')
xlabel('time [s]')


elseif(figure==1 && n_channel==1)
figure;
subplot(4,3,1)
plot(t_vec,C3);
title('C3')
xlabel('time [s]')
subplot(4,3,2)
plot(t_vec,C4);
title('C4')
xlabel('time [s]')
subplot(4,3,3)
plot(t_vec,Cz);
title('Cz')
xlabel('time [s]')
subplot(4,3,4)
plot(t_vec,Fc5);
title('Fc5')
xlabel('time [s]')
subplot(4,3,5)
plot(t_vec,Fc1);
title('Fc1')
xlabel('time [s]')
subplot(4,3,6)
plot(t_vec,Fc2);
title('Fc2')
xlabel('time [s]')
subplot(4,3,7)
plot(t_vec,Fc6);
title('Fc6')
xlabel('time [s]')
subplot(4,3,8)
plot(t_vec,Cp6);
title('Cp6')
xlabel('time [s]')
subplot(4,3,9)
plot(t_vec,Cp2);
title('Cp2')
xlabel('time [s]')
subplot(4,3,10)
plot(t_vec,Cp1);
title('Cp1')
xlabel('time [s]')
subplot(4,3,11)
plot(t_vec,Cp5);
title('Cp5')
xlabel('time [s]')
end

end
