function [data_filtered] = CSP(myData)
% Laplacian performs a spatial filtering following the Common Spatial Pattern method

% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
%
% In the param structure, the sample frequency is saved
% 
% The CSP function returns the filtered EEG signals

load('CSPmatrix.mat');

data_filtered=myData; % So that data_filtered is a structure

% Reorganize as C3 Cz C4
tempTrial=myData.trial{1}(2,:);
myData.trial{1}(2,:)=myData.trial{1}(3,:);
myData.trial{1}(3,:)=tempTrial;

data=myData.trial{1};

multData=W_Csp'*data;

data_filtered.trial{1}=multData;
end
