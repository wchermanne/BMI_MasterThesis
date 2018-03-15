function [outputArg1] = FtLDA(trainedModel,featureIn)
%This function predicts the class of a given feature vector based on a
%linear discriminant analysis.
% Input :The Trained model obtained using the Classification Learner app
% from MATLAB trainedModel; The feature vector featureIn;
%
% Output : The class which the feature vector belongs to 
%
% Example : [MyPrediction] = FtLDA(trainedModel,featureIn);
if(size(featureIn,2) ==1)
elseif(size(featureIn,1) ==1)
    featureIn = (featureIn.');
else
    error('The feature vector is not a column/row vector but a matrix');
end
classTag = trainedModel.ClassificationDiscriminant.ClassNames; 
numOfClass = length(classTag);
if(numOfClass == length(featureIn))
else
    error('The feature vector must have the same format as the trained model');
end

Pk = trainedModel.ClassificationDiscriminant.Prior; %% A priori probabiliy of a class
SigmaK = trainedModel.ClassificationDiscriminant.Sigma; %% Inter-Class cross correlation matrix 
muK = (trainedModel.ClassificationDiscriminant.Mu).'; %% mean of a given class; each row is for a class and each column is for a predictor 

Pxk = ones(numOfClass,1); %% (P(X|K) probability of X given its class. See Bayes 
Pkx = ones(numOfClass,1); 
Px = 0;
for l=1:1:numOfClass
    Pxk(l) = 1/sqrt((2*pi*det(SigmaK))) * exp(-0.5*(featureIn - muK(:,l))*(inv(SigmaK))*(featureIn - muK(:,l)));
    Px = Px + Pxk(l)*Pk(l);
    
    Pkx(l) = Pxk(l)*Pk(l);
end
Pkx = Pkx./(Px);

cost = trainedModel.ClassificationDiscriminant.Cost; %% cost matrix. By default we'll use the complementary matrix of the identity 
cost_tot = cost*Pkx;
[minValue, minIndx] = min(cost_tot);
outputArg1 = classTag(minIndx);
end

