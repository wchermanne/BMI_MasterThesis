function [data_filtered] = CAR(myData)
%% Informations
% This functions applies a Common Average Reference filter to the data
%
% INPUTS 
%
% myData is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% OUTPUTS
%
% data_filtered is the FieldTrip Structure containing the filtered data

%% Code
data=myData.trial{1};
data_filtered=myData;

C3=data(1,:);
C4=data(2,:);
Cz=data(3,:);

C3_prime = C3- (1/2)*(C4+Cz) ;
C4_prime = C4- (1/2)*(C3+Cz) ;
Cz_prime = Cz- (1/2)*(C4+C3) ;


data_filtered.trial{1}=[C3_prime;C4_prime;Cz_prime]; % Returns the ft_structure

end



