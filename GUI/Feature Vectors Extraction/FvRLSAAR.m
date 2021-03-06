function [outputArg1,outputArg2,outputArg3] = FvRLSAAR(y,X,E,theta,Pt,Fs,p,q,lambda)
%This function returns a feature vector using the AAR method based on the RLS algorithm. The current
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
% Example : [ARValue] =FvAutoRegressive(X,t,256,10);
if (lambda <=0 || lambda >=1 )
    error('Choose a valid lambda value for the RLS algorithm (0 < lambda < 1)')
end


if (Fs <=0 || Fs >= 2001)
    error('Choose a valid sampling frequency (1 -> 2KHz)')
end

if(p == size(X,2))
    
elseif(p == size(X,1))
    X = (X.'); %% transpose for having the required format
    if(p == size(X,2))
    else
        error('X does not have the required format');
    end
else
    error('X must has the same length as the order p of the ARMA model');
end

if(q == size(E,2))
    
elseif(q == size(E,1))
    E = (E.'); %% transpose for having the required format
    if(q == size(E,2))
    else
        error('E does not have the required format');
    end
else
    error('E must has the same length as the order q of the ARMA model');
end

if(size(theta,2) == (p+q))
    
elseif(size(theta,1) == (p+q))
    theta = (theta.'); %% transpose for having the required format
    if(size(theta,2) == (p+q))
    else
        error('Theta does not have the required format');
    end
else
    error('Theta must has the same length as the order q+p of the ARMA model');
end

if (size(Pt,1) ~= (p+q) || size(Pt,2) ~= (p+q))
    error('Pt does not have the correct format. It must be a square matrix of size p+q')
end

% RLS algorihtm

Xt = [X,E].';
e = y - theta*Xt;
denom_factor = (Xt.')*Pt*Xt;
Kt = Pt*(Xt./(lambda+denom_factor));

theta = theta + e*Kt;
Pt = (1/lambda)*(eye(p+q)-Kt*(Xt.'))*Pt;


outputArg1 = theta;
outputArg2 = Pt;
outputArg3 = e;
end







