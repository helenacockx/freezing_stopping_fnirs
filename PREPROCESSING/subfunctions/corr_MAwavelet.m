function [data_corr_rs] = corr_MAwavelet(data)
% CORR_MAWAVELET
% This function corrects the motion artifacts with a wavelet function. This
% is an adapted version of the Homer 3 function:
% hmrR_MotionCorrectWavelet. It first downsamples the data to 10 Hz to
% fasten the processing and then reupsamples it to 60 Hz.
%
% INPUT:
% - data: the fNIRS data in fieldtrip structure
%
% OUTPUT:
% - data_corr_rs: the fNIRS data after the wavelet correction
%
% DEPENDENCIES
% - Homer3 (to load db2 wavelet)
%
%%
% chose parameters
iqr = 0.8; % cfr. di Lorenzo et al, 2019
vis = false;

% resample to make faster
cfg = [];
cfg.resamplefs = 10;
data_rs = ft_resampledata(cfg, data);

% pad the data to avoid edge artifacts
data_pad = ft_preproc_padding(data_rs.trial{1}, 'mirror', 60*data_rs.fsample);

dod         = data_pad';
dodWavelet  = dod;

SignalLength = size(dod,1); % #time points of original signal
N = ceil(log2(SignalLength)); % #of levels for the wavelet decomposition
DataPadded = zeros(2^N,1); % data length should be power of 2

% Must provide getAppDir function which
db2path = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\Homer3\Install\db2.mat';

fprintf('Loading %s\n', db2path);
load(db2path);  % Load a wavelet (db2 in this case)

qmfilter = qmf(db2,4); % Quadrature mirror filter used for analysis
L = 4;  % Lowest wavelet scale used in the analysis
for ii = 1:length(data.label)
  DataPadded(1:SignalLength) = dod(:,ii);  % zeros pad data to have length of power of 2
  DataPadded(SignalLength+1:end) = 0;
  
  DCVal = mean(DataPadded);
  DataPadded = DataPadded-DCVal;    % removing mean value
  
  [yn, NormCoef] = NormalizationNoise(DataPadded',qmfilter);
  
  StatWT = WT_inv(yn,L,N,'db2'); % discrete wavelet transform shift invariant
  
  ARSignal = WaveletAnalysis(StatWT,L,'db2',iqr,SignalLength);  % Apply artifact removal
  ARSignal = ARSignal/NormCoef+DCVal;
  
  dodWavelet(:,ii) = ARSignal(1:length(dod));
  
  if vis
    figure; hold on
    plot(dod(:,ii));
    a = ylim();
    plot(dodWavelet(:,ii))
    ylim([a(1)-0.1 a(2)+0.1])
    title(data_rs.label{ii})
  end
end

data_corr = data_rs;
data_corr.trial{1} = ft_preproc_padding(dodWavelet', 'remove', 60*data_rs.fsample); % remove the padding

% remove padding
cfg = [];
data_corr = ft_preprocessing(cfg, data_corr);

% upsample to again match the motion data
cfg = [];
cfg.time = data.time;
data_corr_rs = ft_resampledata(cfg, data_corr);
data_corr_rs.fsample = round(data_corr_rs.fsample);
