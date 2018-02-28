function [data_filtered] = ICA(myData)
% Laplacian performs a spatial filtering following the Independent Component Analysis method

% The input data should be organised in a nxm matrix where n is the number
% of channels and m is the number of samples
%
% In the param structure, the sample frequency is saved
% 
% The ICA function returns the filtered EEG signals
channels=myData.trial{1};
time=myData.time{1};
cfg=[];
cfg.method='runica';
comp=ft_componentanalysis(cfg,myData);
assignin('base', 'comp', comp);

myData.trial{1}=comp.trial{1};
data_filtered=myData;

    %% Independent Comp Analysis (Back to spatial filtering) with Matthieu's method
    ICASig = fastica([myData.trial{1}(1,:); myData.trial{1}(3,:); myData.trial{1}(2,:)], 'numOfIC', 3);
    S=ICASig;
    comp=ICASig;
    
    figure
    subplot(3,1,1)
    plot(time,S(1,:))
    title('S1 Channel with Fast ICA')
    subplot(3,1,2)
    plot(time,S(2,:))
    title('S2 Channel with Fast ICA')
     subplot(3,1,3)
    plot(time,S(3,:))
    title('S3 Channel with Fast ICA')

end
