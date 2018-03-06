function [FeatureVector] = myFeatureExtraction(myData,method)
%% Informations
% This functions builds the feature vector of the data contained in the
% structure myData
%
% INPUTS 
%
% myData is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% method is a string containing the  name of the method to use.
%
% OUTPUTS
%
% FeatureVector is the vector containing the features

%% Code
if strcmp(method,'PeakDetection')
    [tempVector,time]=FvPeakDetection(myData.trial{1},myData.time{1},'max');
elseif strcmp(method,'PowerBands')
    %% To complete
elseif strcmp(method,'Variance')
    [tempVector] = FvVariance(myData.trial{1},myData.time{1},myData.fsample);
elseif strcmp(method,'Wavelets')    
    %% To complete
elseif strcmp(method,'LogNormal')
    %% To complete
elseif strcmp(method,'TemplateMatching')
    [tempVector,time] = FvTemplateMatching(myData.trial{1},myData.time{1},template,myData.fsample)
elseif strcmp(method,'TimeIntegration')
    [tempVector] = FvTimeIntegration(myData.trial{1},myData.time{1},myData.fsample)
elseif strcmp(method,'AR')
    [tempVector]=FvAutoRegressive(myData.trial{1},myData.time{1},myData.fsample,10)
end

FeatureVector=tempVector;
end
