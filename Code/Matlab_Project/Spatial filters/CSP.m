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

data=myData.trial{1};

multData=W_Csp'*data;

data_filtered.trial{1}=multData;
end
