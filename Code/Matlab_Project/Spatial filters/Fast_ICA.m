function [S] = Fast_ICA(t,X,NumOfIC,NumberOfChannel,Nonlinearity,NL)
%The format of X = R^{Ch x N}. This is a short fat matrix where every
%row corresponds to a channel. The algortihm implements a FAST ICA. The
%output is S = R^{S_j * N} (short fat matrix). Notice that the correct way
%for calling this function is e.g. :
% --------------------------------
% S = Fast_ICA(t,X,'NumOfIC',2,'Nonlinearity','tanh')

% Check the number of channel required
if(NumberOfChannel <= 1)
    error('The Number of channel must be > 1');
end

if(length(t) == size(X,2))
    
elseif(length(t) == size(X,1))
    X = (X.');%% transpose for having the required format
    if(length(t) == size(X,2))
    else
        error('X does not have the required format');
    end
else
    error('Time vector must has the same length as the number of rows of X');
end

% Check the number of channel required
if(NumberOfChannel > size(X,1))
    error('The Number of channel must be less than the original number of EEG signals ');
end

% Pre-whitening part 

mean = sum(X,2)./(size(X,2));
for k =1:1:size(X,1)
    X(k,:) = X(k,:) -mean(k);
end

R_xx = X*(X.');
[E,D] = eig(R_xx);
P = (inv(D))^(1/2);
X_white = E*P*(E.');
X = X_white;
% ICA part 
w_p = rand(size(X,1),1);


S = X;
end

