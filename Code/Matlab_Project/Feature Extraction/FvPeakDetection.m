function [outputArg1] = FvPeakDetection(X,t,type)

%This function returns a feature vector using the peak detection method
% Input : vector or matrix of samples X; time vector
% of the sample time; type refers as the chosen method for peak picking. It
% is either 'min' or 'max'
%
% Output : feature vector made of the max or min value for the current
% block or matrix and its time occurence. Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [peak occurenceTime] = FvPeakDetection(X,256,time,'max');
if strcmp(type,'max')
    op = 1;
elseif strcmp(type,'min')
    op = 2;
else
    error('Choose a valid type : min or max');
end 

if(length(t) == size(X,2))
    
elseif(length(t) == size(X,1))
    X = (X.'); %% transpose for having the required format
    if(length(t) == size(X,2))
    else
        error('X does not have the required format');
    end
else
    error('Time vector must has the same length as the number of rows of X');
end

feature_out = ones(size(X,1),1);
time_out = ones(size(X,1),1);
for k = 1:1:size(X,1)
    if (op == 1)
    [feature_out(k),indx] = max(X(k,:));
    time_out(k) = t(indx);
    else
    [feature_out(k),indx] = min(X(k,:)); 
    time_out(k) = t(indx);
    end
end
outputArg1 = feature_out;
end

