function [] = FeatureClassification(FeatureVector,choices)
%% Informations
% This function performs the classification of a FeatureVector using a
% specified method
%
% INPUTS 
%
% FeatureVector is the vector of features
%

% OUTPUTS
%
%% Code
    FeatureExtractionMethod=choices.FeatureExtractionMethod;
    trainedModel=evalin('base','trainedModel');
    fit = trainedModel.predictFcn(FeatureVector);
    assignin('base','fit',fit);

%% Graphes and prints
    nbOfRight=0;
    nbOfLeft=0;
    nbOfRest=0;
    sequenceRight=0;
    sequenceLeft=0;
    nbOfFrames=length(fit);
    duration=choices.duration;
    timestep=(choices.WindowLength)/1000;
    LeftMovement=0;
    RightMovement=0;
    for i=1:length(fit)
        if(sequenceRight>=3)
            RightMovement=1;
            frame=i;
        elseif (sequenceLeft>=3)
            LeftMovement=1;
            frame=1;
        end
        
        if (fit(i)==1)
            nbOfRight=nbOfRight+1;
            sequenceRight=sequenceRight+1;
            sequenceLeft=0;
        elseif(fit(i)==2)
            nbOfLeft=nbOfLeft+1;
            sequenceLeft=sequenceLeft+1;
            sequenceRight=0;
        elseif(fit(i)==3)
            nbOfRest=nbOfRest+1;
            sequenceLeft=0;
            sequenceRight=0;
        end
    end

    
    
    clc
    fprintf('There are %d frames in total \n',length(fit))
    fprintf('%d frames correspond to left movement \n',nbOfLeft)
    fprintf('%d frames correspond to right movement \n',nbOfRight)
    fprintf('%d frames correspond to rest movement \n',nbOfRest)
    myResult=' ';
    if(LeftMovement==1)
            myResult='Left movement detected!';
            fprintf('Left movement detected!\n')
            
    elseif(RightMovement==1)
            myResult='Right movement detected!';
            fprintf('Right movement detected!\n')
    else myResult='No movement detected!'; 
    end
    assignin('base','myResult',myResult);
end