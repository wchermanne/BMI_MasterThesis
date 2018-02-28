function [data_filtered] = LAP(myData)
%% Informations
% This functions applies a Laplacian filter to the data
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
    
d_C3 = [4 8]; % = [distance_to_C4 distance_to_Cz]
d_C4 = [4 8]; % = [distance_to_C3 distance_to_Cz]
d_Cz = [4 8]; % = [distance_to_C3 distance_to_C4]

% Only for one electrode at a time, so for C3 here
    
g_C3 = (1./d_C3)./(sum(1./d_C3));
g_C4 = (1./d_C4)./(sum(1./d_C4));
g_Cz = (1./d_Cz)./(sum(1./d_Cz));

C3_LAP = C3 - g_C3(1).*C4 - g_C3(2).*Cz;
C4_LAP = C4 - g_C4(1).*C3 - g_C4(2).*Cz; 
Cz_LAP = Cz - g_Cz(1).*C3 - g_Cz(2).*C4; 
    
data_filtered.trial{1}=[C3_LAP;C4_LAP;Cz_LAP]; % Returns the ft_structure

end
