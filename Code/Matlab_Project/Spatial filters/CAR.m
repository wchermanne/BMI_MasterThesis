function [data_filtered] = CAR(data)
% CAR performs a spatial filtering following the Common Average Reference
% Method
% Use as
%   [freq] = ft_freqanalysis(cfg, data)
% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
% In the param structure, the sample frequency is saved
% 
% The CAR functions returns the filtered EEG signals
    C3=data(1,:);
    C4=data(2,:);
    Cz=data(3,:);

    C3_prime = C3- (1/2)*(C4+Cz) ;
    C4_prime = C4- (1/2)*(C3+Cz) ;
    Cz_prime = Cz- (1/2)*(C4+C3) ;
    
    data_filtered=[C3_prime;C4_prime;Cz_prime];
end



