%% CREATE_DATAFRAME_NGE
% Converts the trial data of all normal gait events into a table which can
% be read in by R. It contains the fNIRS data averaged over 4 time windows 
% (baseline, pre, post, and post2) for each channel, trial and participant with
% extra info about each trial.
% (baseline = [-10 -5]; pre = [-3 0]; post = [0 3]; post2 = [7 10])
%
% INPUT:
% - *trails*.mat files: nirs data (*trialsnirs_z.mat) 
% segmented in trials based on the events]
%
% OUTPUT:
% data_nrlgait.csv:  csv-file containing the fNIRS data of all normal gait events
% averaged over 3 time windows (baseline, pre and post) for each channel, 
% trial and participant with extra info about the trial

%% Initialisation
root_dir='C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\';

ID={'PD06', 'PD10', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
  'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
  'PD76', 'PD77', 'PD88', 'PD90', 'PD96',...
  'HC02', 'HC04', 'HC06', 'HC13', 'HC19', 'HC21', 'HC28', 'HC29', 'HC33', 'HC34', 'HC35'...
  'HC41', 'HC42', 'HC57', 'HC60', 'HC66', 'HC67', 'HC68', 'HC76', 'HC90', 'HC91'};
ID=sort(ID);

% create template table
varnames = {'HbO', 'HbR', 'ID', 'group', 'trial', 'trigger', 'walkingvsstanding', 'direction', 'channel', 'ROI', 'hemisphere', 'time'};
vartypes = ['doublenan', 'doublenan', 'cellstr', 'cellstr', 'singlenan', repmat({'cellstr'},1,7)]; 

size_dataframe=table(ID', nan(length(ID),1), 'VariableNames', {'ID', 'size'});
data_up = []; %(unpaired)
tot_trials = 1;

%% loop over subjects
for s=1:length(ID)
  fprintf('====== sub-%s ======= \n', ID{s})
  load(fullfile(root_dir, 'processed', 'final', ['sub-' ID{s}], sprintf('sub-%s_trialsnirs_z.mat', ID{s})));
  sf = trials_nirs.fsample;
  
  % only select trials containing normal gait events
  cfg = [];
  cfg.trials = trials_nirs.trialinfo(:,1)==2;
  trials_nge = ft_selectdata(cfg, trials_nirs); %(nge = normal gait events)
  
  % create empty dataframe
  n_tr = length(trials_nge.trial);
  n_chan = sum(~contains(trials_nge.label, {'a ', 'b ', 'c ', 'd '}))/2;
  sz = n_tr * n_chan * 4; % aftwards remove rows containing nans of bad channels
  data_up_pid = table('Size', [sz length(varnames)], 'VariableNames', varnames, 'Variabletypes', vartypes); % unpaired dataset
  ii = 1; % idx table row
  
  %% loop over trials
  for t=1:n_tr
    [trigger, WvS, direction] = value_lookup(trials_nge.trialinfo(t,2));
    % average over time period 
    time_baseline = [find(trials_nge.time{t}==-10):find(trials_nge.time{t}==-5)];
    time_pre = [find(trials_nge.time{t}==-3):find(trials_nge.time{t}==0)];
    time_post = [find(trials_nge.time{t}==0):find(trials_nge.time{t}==3)];
    time_post2 = [find(trials_nge.time{t}==7):find(trials_nge.time{t}==10)];

    data_baseline = nanmean(trials_nge.trial{t}(:, time_baseline),2);
    data_pre = nanmean(trials_nge.trial{t}(:, time_pre),2);
    data_post = nanmean(trials_nge.trial{t}(:, time_post),2);
    data_post2 = nanmean(trials_nge.trial{t}(:, time_post2),2);
    
    %% loop over channels
    bookkeep_chan = zeros(length(trials_nge.label)*2,1);
    for c=1:length(trials_nge.label)
      % check if the channel was not yet processed & if is not a short
      % channel
      if bookkeep_chan(c)==1 | contains(trials_nge.label{c}, {'a ', 'b ', 'c ', 'd '})
        continue
      end
      
      % find corresponding oxy & deoxy channel
      name_split = strsplit(trials_nge.label{c}, ' ');
      idx_chan = find(startsWith(trials_nge.label, name_split{1}));
      if length(idx_chan) ~= 2
        error
      end
      
      % fill in table
      data_up_pid.ID(ii:ii+3) = ID(s); 
      data_up_pid.group(ii:ii+3) = {ID{s}(1:2)};
      data_up_pid.trial(ii:ii+3) = tot_trials; % for each new trial (also in new patient), use new trial number
      data_up_pid.trigger(ii:ii+3) = {trigger};
      data_up_pid.walkingvsstanding(ii:ii+3) = {WvS};
      data_up_pid.direction(ii:ii+3) = {direction};
      data_up_pid.channel(ii:ii+3) = name_split(1);
      [ROI, hemisphere] = ROI_lookup(name_split(1));
      data_up_pid.ROI(ii:ii+3) = {ROI};
      data_up_pid.hemisphere(ii:ii+3) = {hemisphere};
      
      if endsWith(trials_nge.label{idx_chan(1)}, '[O2Hb]') & endsWith(trials_nge.label{idx_chan(2)}, '[HHb]')
        data_up_pid.HbO(ii) = data_baseline(idx_chan(1));
        data_up_pid.HbR(ii) = data_baseline(idx_chan(2));
        data_up_pid.time(ii) = {'baseline'};
        data_up_pid.HbO(ii+1) = data_pre(idx_chan(1));
        data_up_pid.HbR(ii+1) = data_pre(idx_chan(2));
        data_up_pid.time(ii+1) = {'pre'};
        data_up_pid.HbO(ii+2) = data_post(idx_chan(1));
        data_up_pid.HbR(ii+2) = data_post(idx_chan(2));
        data_up_pid.time(ii+2) = {'post'};
        data_up_pid.HbO(ii+3) = data_post2(idx_chan(1));
        data_up_pid.HbR(ii+3) = data_post2(idx_chan(2));
        data_up_pid.time(ii+3) = {'post2'};
      else
        error
      end   
      
      % update bookkeeping
      bookkeep_chan(idx_chan) = 1;
      ii = ii+4; 
      
    end % loop over channels
    tot_trials = tot_trials + 1;
  end % loop over trials
  
  % remove nan containing rows
  data_up_pid = rmmissing(data_up_pid);
  
  % save
%   save(fullfile(root_dir, 'processed', 'final', ['sub-' ID{s}], 'data_nrlgait_pid.mat'), 'data_up_pid')
  data_up = [data_up; data_up_pid];
  size_dataframe.size(s) = height(data_up_pid);
  
end % loop over subjects

writetable(data_up, fullfile(root_dir, 'processed', 'final', 'data_nrlgait.csv'), 'FileType', 'text') 

%% SUBFUNCTIONS
function [trigger, WvS, direction] = value_lookup(trialinfo)

if trialinfo>100
  trigger = 'turn';
  turninfo = num2str(trialinfo);
  if strcmp(turninfo(2), '1')
    direction = 'right';
  elseif strcmp(turninfo(2),'2')
    direction = 'left';
  end
  if strcmp(turninfo(3), '1')
    WvS = 'walking';
  elseif strcmp(turninfo(3), '2')
    WvS = 'standing';
  end
elseif trialinfo>20
  trigger = 'door';
  direction = 'NA'; 
  if trialinfo == 21
    WvS = 'walking';
  elseif trialinfo == 22
    WvS = 'standing';
  end
elseif trialinfo >= 3 & trialinfo <4
  trigger = 'stop';
  if trialinfo == 3.1
    direction = 'turn';
  elseif trialinfo == 3.2
    direction = 'door';
  else
    direction = 'NA';
  end
  WvS = 'walking';
elseif trialinfo == 4.1 | trialinfo == 4
  trigger = 'start';
  direction = 'NA';
  WvS = 'standing';
elseif trialinfo == 4.2
  trigger = 'start_signal';
  direction = 'NA';
  WvS = 'standing';
end

end

function [ROI, hemisphere] = ROI_lookup(channel)
if contains(channel, {'Rx1-Tx2', 'Rx1-Tx3', 'Rx3-Tx2'})
  ROI = 'M1';
  hemisphere = 'right';
elseif contains(channel, {'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7'})
  ROI = 'M1';
  hemisphere = 'left';
elseif contains(channel, {'Rx3-Tx3', 'Rx3-Tx5', 'Rx4-Tx5'})
  ROI = 'PMC';
  hemisphere = 'right';
elseif contains(channel, {'Rx7-Tx8', 'Rx7-Tx10', 'Rx8-Tx10'})
  ROI = 'PMC';
  hemisphere = 'left';
elseif contains(channel, {'Rx4-Tx3', 'Rx4-Tx4', 'Rx2-Tx4'})
  ROI = 'SMA';
  hemisphere = 'right';
elseif contains(channel, {'Rx8-Tx8', 'Rx8-Tx9', 'Rx6-Tx9'})
  ROI = 'SMA';
  hemisphere = 'left';
elseif contains(channel, {'Rx13-Tx17', 'Rx13-Tx18', 'Rx15-Tx17', 'Rx15-Tx18'})
  ROI = 'dlPFC';
  hemisphere = 'right';
elseif contains(channel, {'Rx9-Tx12', 'Rx9-Tx13', 'Rx11-Tx13', 'Rx11-Tx12'})
  ROI = 'dlPFC';
  hemisphere = 'left';
elseif contains(channel,  {'Rx14-Tx19', 'Rx16-Tx19', 'Rx16-Tx20'})
  ROI = 'PPC';
  hemisphere = 'right';
elseif contains(channel,  {'Rx12-Tx14', 'Rx12-Tx15', 'Rx10-Tx14'})
  ROI = 'PPC';
  hemisphere = 'left';
end

end