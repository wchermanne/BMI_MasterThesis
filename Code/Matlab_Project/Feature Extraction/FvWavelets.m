function [outputArg1] = FvWavelets(X,t,Fs,p,type)
%This function returns a feature vector using the Wavelet method.
% Input : vector or matrix of samples X;  time vector
% of the sample time; The sampling frequency Fs, the order of the wavelet model
% p;
%
% Output : feature vector made of the correlation value for the current
% block or matrix with a given template . Notice that if the input X is a
% N*T matrix, the feature vector will have a size of N.
%
% Example : [WaveletValue] =FvWavelets(X,t,Fs,p,type)

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

% wavelet decomposition analysis
wname = type;
CMatrix = [];
LevelMatrix = [];
for ch = 1:1:size(X,1)
    [C,Level] = wavedec(X(ch,:),p,wname);
    CMatrix =[CMatrix ; C];
    LevelMatrix = [LevelMatrix ; Level];
end

wavelet_feature_vec = [];
wavelet_var_tot = 0;

CMatrix = CMatrix(:,LevelMatrix(1,1)+1:end); %% Remove the A_level (Approximation level) component that does not interest us

for kl = 2:1:p+1
    current_length = LevelMatrix(1,kl);
    Cdec = CMatrix(:,1:current_length);
    CMatrix =  CMatrix(:,current_length+1:end); 
    wavelet_tot_current = 0;
    for chX = 1:1:size(X,1)
        wavelet_tot_current = wavelet_tot_current + var(Cdec(chX,:));
    end
    wavelet_var_tot = wavelet_var_tot + wavelet_tot_current;
    featureVec = var(Cdec.');
    wavelet_feature_vec = [wavelet_feature_vec; featureVec.'];
end

outputArg1 = log(wavelet_feature_vec./wavelet_var_tot);
end

