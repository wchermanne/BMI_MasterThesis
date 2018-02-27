function [data_filtered] = LAP(data)
% Laplacian performs a spatial filtering following the Laplacian method

% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
%
% In the param structure, the sample frequency is saved
% 
% The ALAP function returns the filtered EEG signals
    C3=data(1,:);
    C4=data(2,:);
    Cz=data(3,:);
    
    d_C3 = [4 8]; %% Cz followed by C3, distance in cm
    w_hi = (1./d_C3)./(sum(1./d_C3));
    C3_LAP = C3 - w_hi(1).*Cz - w_hi(2).*C4;
    C4_LAP = C3 - w_hi(1).*Cz - w_hi(2).*C3; % To modify
    Cz_LAP = C3 - w_hi(1).*C4 - w_hi(2).*C3; % To modify
    
    data_filtered=[C3_LAP;C4_LAP;Cz_LAP];

end
