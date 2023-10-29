function [run] = remove_badchannels(run, ID)
% REMOVE_BADCHANNELS
% This function detects the bad channels and saves this info in the info
% field of each run (run.info). Bad channels were defined as channels with
% a signal quality index (Sapia et al.2020) of less than 2 during > 50% of the
% standing still periods. Additionally the channels that were bad in 
% >50% of the runs, and all short channels that were bad in at least one
% run are marked as bad channels.
%
% INPUT:
% - run: matlab structure containing the raw nirs & motion data of each run
% - ID: participant-ID of the run that is currently loaded
%
% OUTPUT:
% - run: same as input, but with info about the bad channels added to the
% info field (see run.info.bad_channels)
%
% DEPENDENCIES:
% - artinis fieldtrip functions
% - table2FT
%%
[~, ftpath]=ft_version;
addpath(fullfile(ftpath, 'external', 'artinis')); % add artinis functions

vis = false;
method = 'SQI'; % use the signal quality index (instead of the scalp-coupling index)

%% loop over the runs
for r=1:length(run)
  if isempty(run(r).data_nirs)
    continue
  end
  
  % select the resting data
  fs = run(r).data_nirs.data_raw.fsample;
  eventFT = table2FT(run(r).events, fs);
  eventFT = ft_filter_event(eventFT, 'type', 'lslevent'); % only filter out lsl events
  end_rest_idx = find(strcmp({eventFT(:).value}, 'start_block')); % standing period ends at the start of the new block
  begin_rest_idx = end_rest_idx-1; % standing period begins 1 event before the start of a new block
  % checks
  if any(~contains({eventFT(begin_rest_idx).value}, {'start_run', 'stop_block'}))
    error
  end
  duration_rest = [eventFT(end_rest_idx).timestamp] - [eventFT(begin_rest_idx).timestamp];
  if round(duration_rest, -1) ~= [60 repmat(30,[1 length(duration_rest)-1])] % first standing still of each run is 60 s, other standing still periods are 30 s.
    error
  end
  
  % divide the resting time in equal 10 sec pieces
  begin_sample = [eventFT(begin_rest_idx).sample]';
  begin_sample = sort([begin_sample; begin_sample(1)+30*fs]);
  begin_sample = sort([begin_sample; begin_sample(:) + 10*fs; begin_sample+20*fs]);
  end_sample = begin_sample + 10*fs;
  cfg=[];
  cfg.trl = [begin_sample end_sample zeros(length(begin_sample),1)];
  data_rest = ft_redefinetrial(cfg, run(r).data_nirs.data_raw);
  
  %% detect the bad channels
  switch method
    case 'SCI'
      % select good & bad channels (can also not handle nans)
      cfg                 = [];
      cfg.keepchannel   = 'nan';
      [data_good, good_chan, bad_chan] = ft_nirs_scalpcouplingindex(cfg,data_rest);
      
    case 'SQI'
      % with SQI (Sappa et al.)
      cfg=[];
      cfg.channel = {'Rx*'};
      datain = ft_selectdata(cfg, data_rest);
      cfg=[];
      cfg.windowlength = 10; % default = 10
      cfg.threshold = 1.5; % threshold is below 2
      dataout = ft_nirs_signalqualityindex(cfg, datain);
      
      % find the bad channels for each epoch
      bad_chan_run = {};
      for t=1:length(dataout.time)
        bad_chan_trial{t} = dataout.label(find(any(isnan(dataout.trial{t}(:,:)),2)));
        bad_chan_run = [bad_chan_run; bad_chan_trial{t}];
      end
      
      % indicate the channels that were bad in >= 50% of the epoch
      [counts, chan] = groupcounts(bad_chan_run);
      bad_chan = chan(find(counts>=t/2));
      fprintf('Bad channels of run %d: %d long channels and %d short channels \n', r, sum(~contains(bad_chan,{'a ', 'b ', 'c ', 'd '}))/2, sum(contains(bad_chan, {'a ', 'b ', 'c ', 'd '}))/2)
      display(bad_chan);
      
  end
  
  %% visualize data
  if vis
    cfg =[]; % first convert to conc. for visualization
    cfg.target          = {'O2Hb', 'HHb'};
    cfg.channel         = 'nirs';
    data_conc           = ft_nirs_transform_ODs(cfg, data_rest);
    cfg                = [];
    cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
    cfg.viewmode       = 'vertical';
    cfg.continuous = 'no';
    cfg.channel = [1 2];
    cfg.fontsize=5;
    cfg.nirsscale =5;
    cfg.ylim = [-1 1];
    cfg.linecolor = 'rbkk';
    cfg.colorgroups = repmat([1 2],1, length(data_rest.label)/2)+2*ismember(data_rest.label, bad_chan)';
    ft_databrowser(cfg, data_conc);
%     pause();
  end
  
  %% save the info of the bad channels
  bad_chan = cellfun(@strsplit, bad_chan, 'UniformOutput', false);
  bad_chan = cellfun(@(x) x{1}, bad_chan, 'UniformOutput', false);
  run(r).info.bad_channels = unique(bad_chan);
  
end

%% detect overall bad channels
% detects channels that were bad in >50% of the runs, and all
% short channels that were bad
all_bad_chan = {};
empty_runs = 0; 
for r=1:length(run)
  if isempty(run(r).data_nirs)
    empty_runs = empty_runs + 1;
    continue
  end
  all_bad_chan = [all_bad_chan; run(r).info.bad_channels];
end
[counts, chan] = groupcounts(all_bad_chan);
extr_bad_chan = chan(find(counts>=(numel(run)-empty_runs)/2));
if ~isempty(extr_bad_chan)
  warning('These channels were bad in more than half of the runs: %d long channels and %d short channels', sum(~contains(extr_bad_chan,{'a', 'b', 'c', 'd'})), sum(contains(extr_bad_chan, {'a', 'b', 'c', 'd'})))
  display(extr_bad_chan)
  if length(extr_bad_chan)>16
    warndlg('More than half of the long channels were bad.', ID)
  end
end
bad_short_chan = chan(contains(chan, {'a', 'b', 'c', 'd'}) | contains(chan, {'a-', 'b-', 'c-', 'd-'}));
if ~isempty(bad_short_chan)
  warning('These short channels were bad in at least one of the runs:')
  display(bad_short_chan)
  if length(bad_short_chan)>8
    warndlg('More than half of the short channels were bad.', ID)
  end
  if strcmp(ID, 'HC76')
    warning('Only excluding short channels that were bad in more than half of the runs')
    bad_short_chan = extr_bad_chan(contains(extr_bad_chan, {'a', 'b', 'c', 'd'}));
  end
end

% update the info.bad_channels
for r=1:length(run)
  if isempty(run(r).data_nirs)
    continue
  end
  run(r).info.bad_channels = union(run(r).info.bad_channels, [extr_bad_chan; bad_short_chan]);
end


