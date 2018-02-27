function [mySpatialFilteredData]=SpatialFiltering(data,type,band)
signals=data.trial{1};
if strcmp(type,'CAR')
    tempData=CAR(signals);
elseif strcmp(type,'LAP')
    tempData=LAP(signals);
end

data.trial{1}=tempData;
mySpatialFilteredData=data;

time=data.time{1};
channels=mySpatialFilteredData.trial{1};
titleC3='C3 x time spatially filtered for '
myTitleC3=[titleC3 band];
titleC4='C4 x time spatially filtered for '
myTitleC4=[titleC4 band];
titleCz='Cz x time spatially filtered for '
myTitleCz=[titleCz band];
%% Figures
figure;
ax1=subplot(3,1,1)
plot(time,channels(1,:))
xlabel('time [s]')
title(myTitleC3)
ax2=subplot(3,1,2)
plot(time,channels(2,:))
xlabel('time [s]')
title(myTitleC4)
ax3=subplot(3,1,3)
plot(time,channels(3,:))
xlabel('time [s]')
title(myTitleCz)
end
