function [outputArg1,outputArg2] = FvTemplateMatching(X,t,template,Fs)
%This function returns a feature vector using the correlation and template matching method
% Input : vector or matrix of samples X;  time vector
% of the sample time; A template which is a vector or a matrix. Notice that
% the dimension of the template must be the same as the X vector (matrix);
% The samplign frequency Fs
% 
% Output : feature vector made of the correlation value for the current
% block or matrix with a given template . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [XcorrValue] = FvTemplateMatching(X,t,template,256);
%
% Notice that this method assumes a training session 

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

if(size(X,1) == size(template,1))
    
elseif(size(X,1) == size(template,2))
    template = (template.');  %% transpose for having the required format
    if (size(X,1) == size(template,1))
    else
        error('Template does not have the required format')
    end
else
    error('Template must have the same format as X')
end

feature_out = ones(size(X,1),1);
time_out = ones(size(X,1),1);
for k = 1:1:size(X,1)
    [feature_out(k),indx] = max(xcorr(X(k,:),template(k,:)));
    time_out(k) = t(indx);
end
outputArg1 = feature_out;
outputArg2 = ime_out;
end




