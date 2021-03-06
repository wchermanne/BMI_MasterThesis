function [outputArg1] = FvAutoRegressive(X,t,Fs,p)

%This function returns a feature vector using the AR method. The current
% sample might be predicted as a linear combination of
% the previous samples + a prediction error. The algorithm used here is based on Levinson Durbin algorithm
% Input : vector or matrix of samples X;  time vector
% of the sample time; The sampling frequency Fs, the order of the AR model
% p;
%
% Output : feature vector made of the correlation value for the current
% block or matrix with a given template . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [ARValue] =FvAutoRegressive(X,t,256,10);

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

% Autocorrelation part
r=zeros(p+1,size(X,1));
for ch = 1:1:size(X,1)
    for k=1:p+1
        r(k,ch)=sum(X(ch,k:end).*X(ch,1:end-k+1));
    end
end


% a = predicition polynomial coefficients, including the leading 1 (a0)
% Beware for sign: use polynomial coefficients as-is: freqz(1,a(:,frame))
% is transfer function
% k = reflection coefficients
% Eres: variance of the residual

[N,T]=size(r);
a=zeros(N,T);
a(1,:)=1; % Leading 1 (a0)
k=zeros(N-1,T);

for m=1:N-1
    k(m,:)=sum(r(m+1:-1:2,:).*a(1:m,:),1) ./ sum(r.*a,1);
    a(2:m+1,:)=a(2:m+1,:)-k(m*ones(m,1),:).*a(m:-1:1,:);
end
Eres=sum(r.*a,1);

outputArg1 = reshape(a,[],1);
outputArg2 = k;
outputArg3 = Eres;
end






