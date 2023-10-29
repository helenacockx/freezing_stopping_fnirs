function data_out= auxchannel_regression(data_nirs, data_accel, data_gyr, data_orient)
% AUXCHANNEL_REGRESSION
% auxchannelregression performs a OLS regression of the fNIRS data with the long
% channels as the dependent variable and the short channels and movement
% data as regressors. It first z-transforms the data and uses only the 8 first principle 
% components of the short channel data. The residuals of the regression analysis are the output.
%
% INPUT:
% - data_nirs: the fNIRS data in fieldtrip structure (including the long
% and short channels)
% - data_accel: the accelerometer data of the head in fieldtrip structure
% - data_gyr: the gyroscope data of the head in fieldtrip structure
% - data_orient: the orientation data of the head in fieldtrip structure
%
% OUTPUT:
% - data_out: the fNIRS data after the regression analysis (the residuals)
%
% DEPENDENCIES:/
%%
%% separate long from short channels
cfg = [];
cfg.channel = find(contains(data_nirs.label, {'a ', 'b ', 'c ', 'd '})| contains(data_nirs.label, {'a-', 'b-', 'c-', 'd-'}));
data_short = ft_selectdata(cfg, data_nirs);
cfg.channel = find(~contains(data_nirs.label, {'a ', 'b ', 'c ', 'd '}) & ~contains(data_nirs.label, {'a-', 'b-', 'c-', 'd-'}));
data_long = ft_selectdata(cfg, data_nirs);

%% z-score
% short channels
data_short_z = (data_short.trial{1}-nanmean(data_short.trial{1},2))./nanstd(data_short.trial{1}')';

% IMU
data_accel_z = (data_accel.trial{1}-nanmean(data_accel.trial{1},2))./nanstd(data_accel.trial{1}')';
data_gyr_z = (data_gyr.trial{1}-nanmean(data_gyr.trial{1},2))./nanstd(data_gyr.trial{1}')';
data_orient_z = (data_orient.trial{1}-nanmean(data_orient.trial{1},2))./nanstd(data_orient.trial{1}')';

%% perform PCA
% short channels
[coeff,score,latent,tsquared,explained] = pca(data_short_z'); % assuming that there is only 1 trial!
display(sum(explained(1:8)))
display(explained(1:8)')
short_pca = score(:,1:8);

%% OLS solution
y = data_long.trial{1}(:,:)'; % assuming that there is only 1 trial!

AUX = [short_pca data_accel_z' data_gyr_z' data_orient_z'];
AUX = [repmat(1, size(AUX, 1), 1) AUX];
beta = AUX\y;
yhat = AUX*beta;
res  = y - yhat;

%% update data structure
data_out = data_long;
data_out.trial{1} = res';

