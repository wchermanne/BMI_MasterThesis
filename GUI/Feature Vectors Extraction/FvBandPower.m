function [outputArg1] = FvBandPower(X,t,Fs)
%This function returns a feature vector using the time integration method
%based on Darboux sum.
% Input : vector or matrix of samples X;  time vector
% of the sample time; Sampling frequency Fs;
%
% Output : feature vector made of the integration value for the current
% block or matrix . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [IntValue] = FvTimeIntegration(X,t,256);
if (Fs <=0 || Fs >= 2001)
    error('Choose a valid sampling frequency (1 -> 2KHz)')
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
for k = 1:1:size(X,1)
    feature_out(k) = (X(k,:)*ones(size(X,1),1)).*(1/Fs);
end
outputArg1 = feature_out;
end


