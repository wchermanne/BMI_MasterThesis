function [outputArg1,outputArg2,outputArg3] = FvRAR(Y,A,Fs,p,lambda,Ak)
%This function returns a feature vector using the AAR method based on the RAR algorithm. The current
% sample might be predicted as a linear combination of
% the previous samples + a prediction error. Here the AR coefficient change
% over time (foreach new sample)
% The RAR algortihm is the acronym for Recursive AR.
% Input : vector or matrix of samples Y;  time vector
% of the sample time; The sampling frequency Fs, the order of the AAR model
% p; The update coefficient lambda; the signal covariance Ak, The vector
% for AR coefficients A
%
% Output : feature vector made of the correlation value for the current
% block or matrix with a given template . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [outputArg1,outputArg2,outputArg3] = FvRAR(Y,A,Fs,p,lambda,Ak)
if (lambda <=0 || lambda >=200 )
    error('Choose a valid lambda value for the LMS algorithm (0 < lambda < 200)')
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

if (size(Ak(:,:,1)) ~= [p p])
    error('Ak does not have the expected size (p x p)')
end

% Computing the prediction error
yk = Y(:,end);
Y = Y(:,1:end-1); %% Removing the current sample from Y
ek = ones(size(Y,1),1);
Aknew = ones(p,p,size(Y,1));
% Computing the one step ahead prediction error
% Computing Ak based on Ak-1
for indx = 1:1:size(Y,1) %% Looping for each channel
    ek(indx) = yk(indx) - A(indx,:)*Y(indx,:).';
    Aknew(:,:,indx) = (1-lambda).*Ak(:,:,indx) + lambda.*(Y(indx,:).'*Y(indx,:));
end

Kk = ones(size(Y,1),p);
% Computing Kk based on Ak
for indx = 1:1:size(Y,1) %% Looping for each channel
    denum = lambda*(Y(indx,:)*Aknew(:,:,indx)*Y(indx,:).') +1;
    Kk(indx,:) = (lambda*Aknew(:,:,indx)*Y(indx,:).')./denum;
end


% Computing the new estimate of AR coefficients
Anew = ones(size(A));
for indx = 1:1:size(Y,1) %% Looping for each channel
    Anew(indx,:) = A(indx,:) + Kk(indx,:)*ek(indx);
end

outputArg1 = Anew;
outputArg2 = Aknew;
outputArg3 = ek;
end









