function [outputArg1] = FvTimeIntegration(X,t,Fs)
%% Informations
% This function returns a feature vector using the time integration method
% based on Darboux sum.
%
% INPUTS 
%
% X is the vector or matric of samples of interest
%
% t is the time vector
%
% Fs is the sample frequency
%
% OUTPUTS
%
% outputArg1 is the feature vector

%% Code
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
    feature_out(k) = (X(k,:).^2*ones(size(X,2),1)).*(1/Fs);
end

%% Output
outputArg1 = feature_out;
end


