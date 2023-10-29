function [artifact] = detect_MA(data, bad_chan)
% DETECT_MA
% This detects motion artifacts with a moving average and is an adapted version from Homer 3: hrmR_MotionArtefactByChannel
% The parameters are defined at the beginning of the function (tMotion =
% 0.5; tMask = 1; std_thresh = 65; amp_thresh = 0.05)
%
% INPUT:
% - data: the fNIRS data in fieldtrip structure
% - bad_chan: cell array with the bad channels of the given data
%
% OUTPUT:
% - artifact: artifact fieldtrip structure containing all the detected
% motion artifacts
%
%%
vis = 0; % visualize output?

d = data.trial{1}';
fs = data.fsample;

tMotion = 0.5; 
tMask = 1;
std_thresh = 	65;
amp_thresh = 0.05;

% set artifact buffer for tMask seconds on each side of spike
art_buffer = round(tMask*fs); % time in seconds times sample rate

% loop over channels
bookkeep = zeros(length(data.label),1);
artifact = [];
for i = 1:length(data.label)
  % check if the channel was not yet processed
  if bookkeep(i)==1
    continue
  end
  
  % search for the corresponding channel
  name_split = strsplit(data.label{i}, ' ');
  idx_chan = find(startsWith(data.label, name_split{1}));
  if length(idx_chan) ~= 2
    error
  end
  
  % check if this is not a bad channel
  if contains(data.label{i}, bad_chan)
    bookkeep(idx_chan) = 1;
    continue
  end
  
  % calculate std_diff for each channel
  std_diff = std(d(2:end,idx_chan)-d(1:end-1,idx_chan),0,1);
  
  % calculate max_diff across channels for different time delays
  max_diff = zeros(size(d,1)-1,length(idx_chan));
  for ii=1:round(tMotion*fs)
    max_diff=max([abs(d((ii+1):end,idx_chan)-d(1:(end-ii),idx_chan)); zeros(ii-1,length(idx_chan))], max_diff);
  end

  % find indices with motion artifacts based on std_thresh or amp_thresh
  bad_inds = zeros(size(max_diff));
  mc_thresh=std_diff*std_thresh;
  for ii=1:length(idx_chan)
    bad_inds(:,ii) = max( [max_diff(:,ii)>mc_thresh(ii) max_diff(:,ii)>amp_thresh], [],2);
  end
  bad_inds = find(max(bad_inds,[],2)==1);
  
  % Eliminate time points before or after motion artifacts
  if ~isempty(bad_inds)
    bad_inds=repmat(bad_inds, 1, 2*art_buffer+1)+repmat(-art_buffer:art_buffer,length(bad_inds), 1);
    bad_inds=bad_inds((bad_inds>0)&(bad_inds<=(size(d, 1)-1)));
    
    % Save the MA in an artifact structure
    boolean = zeros(1, length(data.time{1}));
    boolean(bad_inds+1) = 1; % bad inds calculated on diff so add 1
    art_begin = find(diff(boolean)==1);
    art_end = find(diff(boolean)==-1);
    if (length(art_begin) == length(art_end)+1) & boolean(end)==1 % if the data ends with an artifact
      art_end = [art_end length(boolean)];
    end
    art_chan = [art_begin' art_end'];
    
    % visualize the MA
    if vis
    cfg                = [];
    cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
    cfg.viewmode       = 'vertical';
    cfg.continuous     = 'yes';
    cfg.blocksize  = 100;
%     cfg.nirsscale =5;
    cfg.channel = idx_chan;
    cfg.artfctdef.zvalue.artifact = art_chan;
    ft_databrowser(cfg, data);
    end
    
    for ii = 1:length(idx_chan) % add the channel idx (! this must be the same as used for later correction)
      artifact = [artifact; art_chan repmat(idx_chan(ii), size(art_chan,1), 1)];
    end
  end
  
  % mark channels as done
  bookkeep(idx_chan) = 1;
  
end % loop over channels

  
