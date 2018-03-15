function [windowedSignal] = windowing_framing(Fs,t,X,winlen_ratio,winoverlap)
%This function Takes a matrix/vector and windows it using a Hamming window.
% initialize the signal time segment index
% winlen & winoverlap are expressed in milliseconds
% X = R^{N*Ch}, tall thin matrix
% winlen_ratio defines the analysis window length as a multiple of the
% overlap window length according to the follwoing relationship : winlen =
% winlen_ratio * winoverlap (assuming that winlen_ratio is an integer)
% Notice that the last samples are discared for having an integer number of
% frames.
if(Fs <= 0 || winlen_ratio <= 0 || winoverlap <= 0)
    error('Fs & winoverlap & winlen_ratio must be stricly positive')
end
hop = floor((winoverlap./1000)*Fs); %% window overlap in samples
if (~mod(winlen_ratio,1) == 1)
else
    error('Winlen_ration must be an integer');
end

if(length(t) == size(X,1))
    
elseif(length(t) == size(X,2))
    X = (X.');%% transpose for having the required format
    if(length(t) == size(X,1))
    else
        error('X does not have the required format');
    end
else
    error('Time vector must has the same length as the number of rows of X');
end
wlen = hop*winlen_ratio; %% window length in samples
win = hamming(wlen, 'periodic');

frame_numb = floor(size(X,1)./hop);
sample_numb_removing= size(X,1)-(frame_numb)*hop;
X = X(1:end-sample_numb_removing,:);
t = t(1:end-sample_numb_removing);

window_numb = frame_numb - winlen_ratio +1;
if(window_numb <=0)
    error('The analysis window size is too large');
end
for l = 1:1:size(X,2)
    Xw_matrix = zeros(wlen,window_numb);
    t_matrix = zeros(wlen,window_numb);
    indx = 0;
    for k = 1:1:window_numb
        % windowing
        Xw = win.*(X(indx+1:indx+wlen,l));
        % Matrix with time segements
        Xw_matrix(:,k) = Xw;
        if(l==1)
            time = t(indx+1:indx+wlen);
            t_matrix(:,k) = time;
        else
        end
        % update the index
        indx = indx + hop;
    end
    windowedSignal.Channel(l).X = Xw_matrix.';
    if(l==1)
        windowedSignal.time = t_matrix.';
    else
    end
end

end

