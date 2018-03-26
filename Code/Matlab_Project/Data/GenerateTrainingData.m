function [] = GenerateTrainingData(config)

%% Go Through all trials of the database
trainModel=[];
assignin('base','trainModel',trainModel);
config=evalin('base','config');

%% We simulated the signals

if(config.LoadSignals==1)
    nbOfTrials=7; % To modify later
elseif(config.SimulateSignals==1)
    nbOfTrials=20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:nbOfTrials
        %% Read file and create structure
        if(config.LoadSignals==1)
        myFilename=['testRight' num2str(i) '.txt'];
        assignin('base','myFilename',myFilename);
        [data,time] = ReadData(myFilename,3,250,3,0);
        if(config.NbOfElectrodes==3)
            rawData=data(1:3,:);
            myTrainingData.label={'C3' 'CZ' 'C4'};
        elseif(config.NbOfElectrodes==11)
            rawData=data(1:11,:);
            myFtData.label={'C3' 'CZ' 'C4' 'Fc5' 'Fc1' 'Fc2' 'Fc6' 'Cp5' 'Cp1' 'Cp2' 'Cp6'};
        end
        myTrainingData=[];
        myTrainingData.fsample=250;
        myTrainingData.trial={rawData};
        myTrainingData.sampleinfo=[1 length(data(1,:))];
        myTrainingData.time={time};
        myTrainingData.cfg=[];

        elseif(config.SimulateSignals==1)
        myTrial=load(['EEG_Signals_Right_Trial_' num2str(i) '.mat']);
        assignin('base', 'myTrial',myTrial)
        %% Create the FieldTrip structure
        fsample=myTrial.Fs;
        rawData=[myTrial.C3; myTrial.Cz; myTrial.C4];
        myTrainingData=[];
        myTrainingData.fsample=fsample;
        myTrainingData.trial={rawData};
        myTrainingData.label={'C3' 'CZ' 'C4'};
        myTrainingData.sampleinfo=[1 length(myTrial.C3)];
        myTrainingData.time={myTrial.time};
        myTrainingData.cfg=[];
        end
        %% Pre-processing
        myTemporalFilteredTrainingData=temporalFtFiltering(myTrainingData,config.LFreq,config.HFreq,config.NotchFilter,config.FilterType,config.winType,config.NFilter,' ',0);

        %% Spatial Filtering
        mySpatialFilteredTrainingData=SpatialFiltering(myTemporalFilteredTrainingData,config.SpatialFilteringMethod,' ',0);

        %% Feature Extraction
        % Feature extraction avec windowing de 3s chacune. Comme le mouvement
        % se passe tjrs à 10 s, c'est la window entre 9 et 12s qui est
        % intéressante pour le mouvement

        %% Attention au choix des frames, c pas top là!
    if(config.LoadSignals==1)
        [FeatureVectorTrainingData] = myFeatureExtraction(mySpatialFilteredTrainingData,config.FeatureExtractionMethod,1,5000,config.AROrder,config.WaveletType,config.WaveletLevel);
        FeatureVectorMovementRightTrainingData(i,:)=FeatureVectorTrainingData(2,:);
        FeatureVectorRestTrainingData(i,:)=FeatureVectorTrainingData(1,:); 
    elseif(config.SimulateSignals==1)
        [FeatureVectorTrainingData] = myFeatureExtraction(mySpatialFilteredTrainingData,config.FeatureExtractionMethod,1,3000,config.AROrder,config.WaveletType,config.WaveletLevel);
        FeatureVectorMovementRightTrainingData(i,:)=FeatureVectorTrainingData(4,:);
        FeatureVectorRestTrainingData(i,:)=FeatureVectorTrainingData(2,:);    
    end
end
    

        classVector=ones(nbOfTrials,1);
       % [FeatureVectorMovementRightTrainingData]=FvLogNormalisation(FeatureVectorMovementRightTrainingData);
        FeatureVectorMovementRightTrainingData=[FeatureVectorMovementRightTrainingData, classVector];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:1:nbOfTrials
        if(config.LoadSignals==1)
        filename=['testLeft' num2str(i) '.txt'];
        [data,time] = ReadData(filename,3,250,3,0);
        if(config.NbOfElectrodes==3)
            rawData=data(1:3,:);
            myTrainingData.label={'C3' 'CZ' 'C4'};
        elseif(config.NbOfElectrodes==11)
            rawData=data(1:11,:);
            myFtData.label={'C3' 'CZ' 'C4' 'Fc5' 'Fc1' 'Fc2' 'Fc6' 'Cp5' 'Cp1' 'Cp2' 'Cp6'};
        end
        myTrainingData=[];
        myTrainingData.fsample=250;
        myTrainingData.trial={rawData};
       
        myTrainingData.sampleinfo=[1 length(data(1,:))];
        myTrainingData.time={time};
        myTrainingData.cfg=[];

        elseif(config.SimulateSignals==1)
        myTrial=load(['EEG_Signals_Left_Trial_' num2str(i) '.mat']);
        %% Create the FieldTrip structure
        fsample=myTrial.Fs;
        rawData=[myTrial.C3; myTrial.Cz; myTrial.C4];
        myTrainingData=[];
        myTrainingData.fsample=fsample;
        myTrainingData.trial={rawData};
        myTrainingData.label={'C3' 'CZ' 'C4'};
        myTrainingData.sampleinfo=[1 length(myTrial.C3)];
        myTrainingData.time={myTrial.time};
        myTrainingData.cfg=[];
        end
        %% Pre-processing
        myTemporalFilteredTrainingData=temporalFtFiltering(myTrainingData,config.LFreq,config.HFreq,config.NotchFilter,config.FilterType,config.winType,config.NFilter,' ',0);

        %% Spatial Filtering
        mySpatialFilteredTrainingData=SpatialFiltering(myTemporalFilteredTrainingData,config.SpatialFilteringMethod,' ',0);

        %% Feature Extraction
        % Feature extraction avec windowing de 3s chacune. Comme le mouvement
        % se passe tjrs à 10 s, c'est la window entre 9 et 12s qui est
        % intéressante pour le mouvement

    if(config.LoadSignals==1)
        [FeatureVectorTrainingData] = myFeatureExtraction(mySpatialFilteredTrainingData,config.FeatureExtractionMethod,1,5000,config.AROrder,config.WaveletType,config.WaveletLevel);
        FeatureVectorMovementLeftTrainingData(i,:)=FeatureVectorTrainingData(2,:);
        FeatureVectorRestTrainingData(nbOfTrials+i,:)=FeatureVectorTrainingData(1,:); 
    elseif(config.SimulateSignals==1)
        [FeatureVectorTrainingData] = myFeatureExtraction(mySpatialFilteredTrainingData,config.FeatureExtractionMethod,1,3000,config.AROrder,config.WaveletType,config.WaveletLevel);
        FeatureVectorMovementLeftTrainingData(i,:)=FeatureVectorTrainingData(4,:);
        FeatureVectorRestTrainingData(nbOfTrials+i,:)=FeatureVectorTrainingData(2,:);    
    end

    end

%% We loaded the signals

    RestVector=ones(nbOfTrials*2,1)*3;
   % [FeatureVectorMovementLeftTrainingData]=FvLogNormalisation(FeatureVectorMovementLeftTrainingData);
    FeatureVectorMovementLeftTrainingData=[FeatureVectorMovementLeftTrainingData, classVector*2];

   % [FeatureVectorRestTrainingData]=FvLogNormalisation(FeatureVectorRestTrainingData);
    FeatureVectorRestTrainingData=[FeatureVectorRestTrainingData,RestVector];

    FinalTrainingData=[FeatureVectorMovementRightTrainingData;FeatureVectorMovementLeftTrainingData;FeatureVectorRestTrainingData];

    %% Save the useful vectors in the workspace
    assignin('base', 'FinalTrainingData', FinalTrainingData);
    assignin('base', 'FeatureVectorMovementRightTrainingData', FeatureVectorMovementRightTrainingData);
    assignin('base', 'FeatureVectorMovementLeftTrainingData', FeatureVectorMovementLeftTrainingData);
    assignin('base', 'FeatureVectorRestTrainingData', FeatureVectorRestTrainingData);



end
