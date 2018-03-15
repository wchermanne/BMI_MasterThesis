function [] = TrainingData(choices)

%% Go Through all trials of the database
nbOfTrials=20;
for i=1:1:nbOfTrials
    myTrial=load(['EEG_Signals_Right_Trial_' num2str(i) '.mat']);
    assignin('base', 'myTrial',myTrial)
    
    %% Create the FieldTrip structure
    fsample=myTrial.Fs;
    rawData=[myTrial.C3; myTrial.C4; myTrial.Cz];
    myTrainingData=[];
    myTrainingData.fsample=fsample;
    myTrainingData.trial={rawData};
    myTrainingData.label={'C3' 'C4' 'CZ'};
    myTrainingData.sampleinfo=[1 length(myTrial.C3)];
    myTrainingData.time={myTrial.time};
    myTrainingData.cfg=[];
    
    %% Pre-processing
    myTemporalFilteredTrainingData=temporalFtFiltering(myTrainingData,choices.LFreq,choices.HFreq,choices.NotchFilter,choices.FilterType,choices.winType,choices.NFilter,' ',0);

    %% Spatial Filtering
    mySpatialFilteredTrainingData=SpatialFiltering(myTemporalFilteredTrainingData,choices.SpatialFilteringMethod,' ',0);
    
    %% Feature Extraction
    % Feature extraction avec windowing de 3s chacune. Comme le mouvement
    % se passe tjrs à 10 s, c'est la window entre 9 et 12s qui est
    % intéressante pour le mouvement
    [FeatureVector] = myFeatureExtraction(mySpatialFilteredTrainingData,choices.FeatureExtractionMethod,1,3000,choices.AROrder,choices.WaveletType,choices.WaveletLevel);
    FeatureVectorMovementRight(i,:)=FeatureVector(4,:);
    FeatureVectorRest(i,:)=FeatureVector(2,:);
end
    classVector=ones(nbOfTrials,1);
    FeatureVectorMovementRight=[FeatureVectorMovementRight, classVector];
    
    assignin('base', 'FeatureVectorMovementRight', FeatureVectorMovementRight);
for i=1:1:20
    myTrial=load(['EEG_Signals_Left_Trial_' num2str(i) '.mat']);
    
    %% Create the FieldTrip structure
    fsample=myTrial.Fs;
    rawData=[myTrial.C3; myTrial.C4; myTrial.Cz];
    myTrainingData=[];
    myTrainingData.fsample=fsample;
    myTrainingData.trial={rawData};
    myTrainingData.label={'C3' 'C4' 'CZ'};
    myTrainingData.sampleinfo=[1 length(myTrial.C3)];
    myTrainingData.time={myTrial.time};
    myTrainingData.cfg=[];
    
    %% Pre-processing
    myTemporalFilteredTrainingData=temporalFtFiltering(myTrainingData,choices.LFreq,choices.HFreq,choices.NotchFilter,choices.FilterType,choices.winType,choices.NFilter,' ',0);

    %% Spatial Filtering
    mySpatialFilteredTrainingData=SpatialFiltering(myTemporalFilteredTrainingData,choices.SpatialFilteringMethod,' ',0);
    
    %% Feature Extraction
    % Feature extraction avec windowing de 3s chacune. Comme le mouvement
    % se passe tjrs à 10 s, c'est la window entre 9 et 12s qui est
    % intéressante pour le mouvement
    
    [FeatureVector] = myFeatureExtraction(mySpatialFilteredTrainingData,choices.FeatureExtractionMethod,1,3000,choices.AROrder,choices.WaveletType,choices.WaveletLevel);
    FeatureVectorMovementLeft(i,:)=FeatureVector(4,:);
    FeatureVectorRest(20+i,:)=FeatureVector(2,:); % Car il y a déjà eu le rest du right donc je dois commencer 20 indices plus loin

end
    RestVector=ones(nbOfTrials*2,1)*3;
    FeatureVectorMovementLeft=[FeatureVectorMovementLeft, classVector*2];
    FeatureVectorRest=[FeatureVectorRest,RestVector];
    FinalTrainingData=[FeatureVectorMovementRight;FeatureVectorMovementLeft;FeatureVectorRest];
    
    %% Save the useful vectors in the workspace
    assignin('base', 'FinalTrainingData', FinalTrainingData);

end
