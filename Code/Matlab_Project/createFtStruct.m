function [myFtData]=createFtStruct()
%% Load the files and data
load('parameters.mat'); % Parameters of eeg_generator
load('rawData.mat'); % Data generated by eeg_generator
rawData=[C3;C4;Cz]


%% Create the FieldTrip structure
fsample=parameters.fsample;
t_end=parameters.duration;

myFtData=[];
myFtData.fsample=fsample;
myFtData.trial={rawData};
myFtData.label={'C3' 'C4' 'CZ'};
myFtData.sampleinfo=[1 length(C3)];
myFtData.time={parameters.time};
myFtData.cfg=[];
end
