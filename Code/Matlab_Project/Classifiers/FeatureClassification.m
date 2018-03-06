function [] = FeatureClassification(FeatureVector,method)
%% Informations
% This function performs the classification of a FeatureVector using a
% specified method
%
% INPUTS 
%
% FeatureVector is the vector of features
%
% method is the method used for the classification. It can be
% 'KNN','LDA','SVM'
%
% OUTPUTS
%
%% Code
    if strcmp(method,'kNN')
        KNN=load('KNN.mat');
        assignin('base', 'KNN', KNN);
        [outputArg1] = FtNaiveKNN(KNN.trainedModel,FeatureVector);
        assignin('base', 'outputArg1', outputArg1);
    end
end