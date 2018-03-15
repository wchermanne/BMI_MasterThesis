function [FeatureVector] = myFeatureExtraction(myData,method,WindowLengthRatio,WindowLengthOverlap,ARorder,WaveletType,WaveletLevel)
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

%% Windowing and framing first

[myWindowedData] = windowing_framing(myData.fsample,myData.time{1},myData.trial{1},WindowLengthRatio,WindowLengthOverlap);
assignin('base', 'myWindowedData', myWindowedData);

nbOfChannels=size(myWindowedData.channels);
nbOfChannels=nbOfChannels(2);
nbOfFrames=size(myWindowedData.channels(1).X);
nbOfFrames=nbOfFrames(1);
myFeatureVectorMatrix=[];

%% Code for the feature extraction
if strcmp(method,'PeakDetection')
    for(i=1:nbOfChannels) % Loop for the channels
        for(j=1:nbOfFrames) % Loop for the different frames 
            tempVector(j,:)=FvPeakDetection(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:),'max');
            assignin('base', 'tempVector', tempVector);
        end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector); % To concatenate the feture vectors for the different channels
    end
elseif strcmp(method,'PowerBands')
    for(i=1:nbOfChannels) % Loop for the channels
        for(j=1:nbOfFrames) % Loop for the different frames 
            tempVector(j,:)=FvPowerBands(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:));
            assignin('base', 'tempVector', tempVector);
        end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector); % To concatenate the feture vectors for the different channels
    end
elseif strcmp(method,'Variance')
        for(i=1:nbOfChannels)
            for(j=1:nbOfFrames)
            tempVector(j,:)=FvVariance(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:),myData.fsample);
            end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector);
        end
        
        
elseif strcmp(method,'Wavelets')    
        for(i=1:nbOfChannels)
            for(j=1:nbOfFrames)
            tempVector(j,:)=FvWavelets(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:),myData.fsample,WaveletLevel,WaveletType);
            end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector);
        end

elseif strcmp(method,'TemplateMatching')
    %% To complete

    
elseif strcmp(method,'TimeIntegration')
        for(i=1:nbOfChannels)
            for(j=1:nbOfFrames)
            tempVector(j,:)=FvTimeIntegration(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:),myData.fsample);
            end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector);
        end

elseif strcmp(method,'AR')
        for(i=1:nbOfChannels)
            for(j=1:nbOfFrames)
            tempVector(j,:)=FvAutoRegressive(myWindowedData.channels(i).X(j,:),myWindowedData.time(j,:),myData.fsample,ARorder);
            end
        myFeatureVectorMatrix=cat(2,myFeatureVectorMatrix,tempVector);
        end    
end

FeatureVector=myFeatureVectorMatrix;
assignin('base', 'FeatureVector', FeatureVector);

end
