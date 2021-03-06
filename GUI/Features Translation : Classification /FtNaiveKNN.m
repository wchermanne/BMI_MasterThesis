function [outputArg1] = FtNaiveKNN(trainedModel,featureIn)
%This function predicts the class of a given feature vector based on a
% k nearest neighboors approach.
% Input :The Trained model obtained using the Classification Learner app
% from MATLAB trainedModel; The feature vector featureIn;
%
% Output : The class which the feature vector belongs to
%
% Example : [MyPrediction] = FtNaiveKNN(trainedModel,featureIn);
%
% Notice here that the algortihm used here is the simplest one.
% A clussering approach could be further investigated
if(size(featureIn,2) ==1)
elseif(size(featureIn,1) ==1)
    featureIn = (featureIn.');
else
    error('The feature vector is not a column/row vector but a matrix');
end
classTag = trainedModel.ClassificationKNN.ClassNames;
numOfClass = length(classTag);
if(numOfClass == length(featureIn))
else
    error('The feature vector must have the same format as the trained model');
end

numOfNeighbors = trainedModel.ClassificationKNN.NumNeighbors;
dist = trainedModel.ClassificationKNN.Distance;
Y = trainedModel.ClassificationKNN.Y;
X = [trainedModel.ClassificationKNN.X.column_1 trainedModel.ClassificationKNN.X.column_2];
outputArg1 = classTag(minIndx);

distVec = ones(size(X,1),1);
switch dist
    case 'cosine'
        for l = 1:1:size(X,1)
            distVec(l) = 1- (X(l,:)*featureIn)/sqrt((featureIn.'*featureIn)*(X(l,:)*X(l,:).'));
        end
    case 'euclidean'
        for l = 1:1:size(X,1)
            distVec(l) = (X(l,:) - featureIn.')*(X(l,:).' - featureIn);
        end
    otherwise
end

[XSorted,Xindx] = sort(distVec,'ascent');
classNN = Y(Xindx(1:numOfNeighbors));
classNNunique = unique(classNN);
[out,outIndx] = max(histc(classNN,classNNunique));
outputArg1 = classNNunique(outIndx);

end


