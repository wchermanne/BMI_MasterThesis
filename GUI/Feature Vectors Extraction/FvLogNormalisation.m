function [outputArg1] = FvLogNormalisation(feature_in)
%This function modfies the feature vector for normalizing it.
% Input : input feature vector feature_in;
%
% Output : feature vector is log-normalized.
%
% Example : [feature_out] = FvLogNormalisation(feature_in);
sum_fv = sum(feature_in);
feature_reg = feature_in./(sum_fv);
if(sum(feature_reg>0)==0) %% Check that each entry in feature_reg is > 0 
    feature_out = log(feature_reg);
else
    error('The vector is totally/partially negative')
end
outputArg1 = feature_out;
end

