function [features] = power_bands(data,time)
% power_bands performs a feature extraction following the power_bands method

% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
%
% In the param structure, the sample frequency is saved
% 
features= trapz(data.^2)
end
