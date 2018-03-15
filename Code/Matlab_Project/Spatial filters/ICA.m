function [data_filtered] = ICA(myData,figure)
%% Informations
% This functions applies a ICA filter to the data
%
% INPUTS 
%
% myData is the FieldTrip Structure containing the trials and the relevant
% data. The signals are contained in data.trial{1}
%
% OUTPUTS
%
% data_filtered is the FieldTrip Structure containing the filtered data

%% Code using FieldTrip
channels=myData.trial{1};
time=myData.time{1};

% cfg=[];
% cfg.method='runica';
% comp=ft_componentanalysis(cfg,myData);
% assignin('base', 'comp', comp);

%% Independent Comp Analysis (Back to spatial filtering) with Matthieu's method

ICASig = fastica([myData.trial{1}(1,:); myData.trial{1}(3,:); myData.trial{1}(2,:)], 'numOfIC', 3); %C3 Cz C4
S=ICASig;
myData.trial{1}=ICASig;
data_filtered=myData;


% 
%     comp=ICASig;
% if(figure==1)
%     figure
%     subplot(3,1,1)
%     plot(time,S(1,:))
%     title('S1 Channel with Fast ICA')
%     subplot(3,1,2)
%     plot(time,S(2,:))
%     title('S2 Channel with Fast ICA')
%      subplot(3,1,3)
%     plot(time,S(3,:))
%     title('S3 Channel with Fast ICA')
% end
end
