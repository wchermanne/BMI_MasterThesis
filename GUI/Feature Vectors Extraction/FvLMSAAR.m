function [outputArg1,outputArg2,outputArg3] = FvLMSAAR(Y,A,Fs,p,lambda,Vk,type)
%This function returns a feature vector using the AAR method based on the LMS algorithm. The current
% sample might be predicted as a linear combination of
% the previous samples + a prediction error. Here the AR coefficient change
% over time (foreach new sample)
% Input : vector or matrix of samples X;  time vector
% of the sample time; The sampling frequency Fs, the order of the AAR model
% p;
%
% Output : feature vector made of the correlation value for the current
% block or matrix with a given template . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [ARValue] =FvLMSAAR(y,X,E,theta,Pt,256,8,0.2);
if (lambda <=0 || lambda >=200 )
    error('Choose a valid lambda value for the LMS algorithm (0 < lambda < 2)')
end


if (Fs <=0 || Fs >= 2001)
    error('Choose a valid sampling frequency (1 -> 2KHz)')
end

if(p+1 == size(Y,2))
    
elseif(p+1 == size(Y,1))
    Y = (Y.'); %% transpose for having the required format
    if(p+1 == size(Y,2))
    else
        error('Y does not have the required format');
    end
else
    error('X must has the same length as the order p of the AR model');
end

if(p == size(A,2))
    
elseif(p == size(A,1))
    A = (A.'); %% transpose for having the required format
    if(p == size(A,2))
    else
        error('A does not have the required format');
    end
else
    error('A must has the same length as the order p of the AR model');
end

% Computing the prediction error
yk = Y(:,end);
Y = Y(:,1:end-1); %% Removing the current sample from Y
ek = ones(size(Y,1),1);

% Computing the one step ahead prediction error
for indx = 1:1:size(Y,1) %% Looping for each channel
    ek(indx) = yk(indx) - A(indx,:)*Y(indx,:).';
end

if(strcmp(type,'LMS1'))
    % Compunting the MSY
    MSY = var(Y.');
    
elseif(strcmp(type,'LMS2'))
    MSY = (1-lambda).*Vk + lambda.*ek.^2;
    
else
    error('Choose between LMS1 & LMS2')
end

% Computing the new estimate of AR coefficients
Anew = ones(size(A));
for indx = 1:1:size(Y,1) %% Looping for each channel
    Anew(indx,:) = A(indx,:) + (lambda/MSY(indx))*ek(indx)*Y(indx,:);
end

outputArg1 = Anew;
outputArg2 = MSY;
outputArg3 = ek;
end







