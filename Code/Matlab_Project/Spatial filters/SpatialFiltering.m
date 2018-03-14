function [mySpatialFilteredData]=SpatialFiltering(data,type,band,figure)

%% Informations
% This functions applies a spatial filter to the input data
%
% INPUTS 
%
% data is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% type is a string containing the  name of the filter to apply. It can be
% 'CAR', 'LAP' or 'ICA'
%
% band is the string containing the information about the frequency band
% being processed
%
%
% OUTPUTS
%
% mySpatialFilteredData is the FieldTrip Structure containing the filtered
% data
%

%% Code

tempData=data; % tempData is a structure, such as data

if strcmp(type,'CAR')
    tempData=CAR(data);
elseif strcmp(type,'LAP')
    tempData=LAP(data);
elseif strcmp(type,'ICA')
    tempData=ICA(data);
elseif strcmp(type,'CSP')
    tempData=CSP(data);
end

mySpatialFilteredData=tempData; % Assign tempData to the output
assignin('base', 'mySpatialFilteredData', mySpatialFilteredData);





%% Figures
if (figure==1)
titleC3='C3 x time spatially filtered for '
myTitleC3=[titleC3 band];
titleC4='C4 x time spatially filtered for '
myTitleC4=[titleC4 band];
titleCz='Cz x time spatially filtered for '
myTitleCz=[titleCz band];

time=tempData.time{1};
channels=tempData.trial{1};

figure;
ax1=subplot(3,1,1)
plot(time,channels(1,:))
xlabel('time [s]')
title(myTitleC3)
ax2=subplot(3,1,2)
plot(time,channels(2,:))
xlabel('time [s]')
title(myTitleC4)

if (strcmp(type,'CSP')==0)
    
ax3=subplot(3,1,3)
plot(time,channels(3,:))
xlabel('time [s]')
title(myTitleCz)
end
end
end
