function [myFtData]=createFtStruct(rawData,config)

%% Informations
% This functions builds the FieldTripeStrucre
%
%
% OUTPUTS
%
% myFtData is the FieldTrip structure containing the data created by the
% function eeg_generator and stored in parameters.mat and rawData.mat

%% Code

%% Create the FieldTrip structure
t_end=config.duration;

myFtData=[];
myFtData.fsample=config.fsample;
myFtData.trial={rawData};
if(config.NbOfElectrodes==3)
    myFtData.label={'C3' 'CZ' 'C4'};
elseif(config.NbOfElectrodes==11)
    myFtData.label={'C3' 'CZ' 'C4' 'Fc5' 'Fc1' 'Fc2' 'Fc6' 'Cp5' 'Cp1' 'Cp2' 'Cp6'};
end
myFtData.sampleinfo=[1 length(rawData(1,:))];
myFtData.time={config.time};
myFtData.cfg=[];
end
