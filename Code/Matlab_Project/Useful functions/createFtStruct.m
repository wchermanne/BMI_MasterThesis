function [myFtData]=createFtStruct()

%% Informations
% This functions builds the FieldTripeStrucre
%
%
% OUTPUTS
%
% myFtData is the FieldTrip structure containing the data created by the
% function eeg_generator and stored in parameters.mat and rawData.mat

%% Code
%% Load the files and data
rawData =evalin('base','rawData');
EEG_parameters =evalin('base','EEG_parameters');

%% Create the FieldTrip structure
fsample=EEG_parameters.fsample;
t_end=EEG_parameters.duration;

myFtData=[];
myFtData.fsample=fsample;
myFtData.trial={rawData};
myFtData.label={'C3' 'C4' 'CZ'};
myFtData.sampleinfo=[1 length(rawData(1,:))];
myFtData.time={EEG_parameters.time};
myFtData.cfg=[];
end
