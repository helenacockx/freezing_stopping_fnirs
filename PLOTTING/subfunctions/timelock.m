function [data_timelock_BL] = timelock(trials, ID, vis, type, condition, baselinecorr)
%% TIMELOCK
% Computes the timelock average for the given participant and the given
% condition. It baseline corrects the data if asked for.
% 
% INPUT: 
% - trials: the epoched fNIRS data of the given participant
% - ID: participant-ID of the given participant
% - vis: visualize the timelock average of this participant? (1/0)
% - type: type of the condition to calculate the timelock average from
% ('FOG' or 'normal')
% - condition: condition of which to calculate the timelock average from
% (e.g., 'turn R walking', 'stop turn', 'doorway')
% - baselinecorr: whether or not to baseline correct the data with a
% baseline from -10 to -5 s (1/0)
%
% OUTPUT:
% - data_timelock_BL: (baseline corrected) timelock average for the given
% condition of the given participant

%% load layout
load('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\scripts\layout.mat');
layout.width = 2*layout.width;
layout.height = 2*layout.height;

%% retrieve the correct condition
switch type
  case 'FOG'
    trlidx = trials.trialinfo(:,1)==1;
  case 'normal'
    trlidx = trials.trialinfo(:,1)==2;
  otherwise
    error('Type is not found')
end

% condition
if contains(condition, 'turn') & ~contains(condition, 'stop')
  trlidx = trlidx & trials.trialinfo(:,2)>100;
  if contains(condition, 'R')
    trlidx = trlidx & trials.trialinfo(:,2)<120;
    if contains(condition, 'walking')
      trlidx = trlidx & trials.trialinfo(:,2) == 111;
    elseif contains(condition, 'standing')
      trlidx = trlidx & trials.trialinfo(:,2) == 112;
    else
      warning('No differentiation is made between walking and standing');
    end
  elseif contains(condition, 'L')
    trlidx = trlidx & trials.trialinfo(:,2)>120;
    if contains(condition, 'walking')
      trlidx = trlidx & trials.trialinfo(:,2) == 121;
    elseif contains(condition, 'standing')
      trlidx = trlidx & trials.trialinfo(:,2) == 122;
    else
      warning('No differentiation is made between walking and standing');
    end
  elseif contains(condition, 'walking')
    trlidx = trlidx & (trials.trialinfo(:,2) == 111|trials.trialinfo(:,2) == 121);
  elseif contains(condition, 'standing')
    trlidx = trlidx & (trials.trialinfo(:,2) == 112|trials.trialinfo(:,2) == 122);
  else
    warning('No differentiation is made between left and right turns or standing and walking')
  end
elseif contains(condition, 'doorway')
  trlidx = trlidx & trials.trialinfo(:,2) > 20;
  if contains(condition, 'walking')
    trlidx = trlidx & trials.trialinfo(:,2)==21;
  elseif contains(condition, 'standing')
    trlidx = trlidx & trials.trialinfo(:,2)==22;
  else
    warning('No differentiation is made between doorway walking and standing');
  end
elseif contains(condition, 'stop')
  trlidx = trlidx & (trials.trialinfo(:,2) >= 3 & trials.trialinfo(:,2) < 4);
  if contains(condition, 'turn')
    trlidx = trlidx & trials.trialinfo(:,2) == 3.1;
  elseif contains(condition, 'door')
    trlidx = trlidx & trials.trialinfo(:,2) == 3.2;
  end
elseif contains(condition, 'start')
  if contains(condition, 'signal')
     trlidx = trlidx & trials.trialinfo(:,2) == 4.2;
  else
    trlidx = trlidx & (trials.trialinfo(:,2) == 4.1 | trials.trialinfo(:,2) == 4);
  end
else
  warning('No further differentiation of type was given')
end


%% timelock analysis 
cfg = [];
cfg.nanmean = 'yes';
cfg.trials = find(trlidx);
if ~isempty(cfg.trials)
  data_timelock = ft_timelockanalysis(cfg, trials);
else
  data_timelock = {};
end

% baseline correction
if baselinecorr
  cfg = [];
  cfg.baseline = [-10 -5];
  if ~isempty(data_timelock)
    data_timelock_BL = ft_timelockbaseline(cfg, data_timelock);
  else
    data_timelock_BL = {};
  end
else
  data_timelock_BL = data_timelock;
end

%% plot
if vis
    if ~isempty(data_timelock_BL)
    % first separate O2Hb and HHb channels
    cfg=[];
    cfg.channel='* [O2Hb]';
    data_O2Hb=ft_selectdata(cfg, data_timelock_BL);
    % and rename labels such that they have the same name as HHb channels
    for i=1:length(data_O2Hb.label)
      tmp = strsplit(data_O2Hb.label{i});
      data_O2Hb.label{i}=tmp{1};
    end
    cfg=[];
    cfg.channel='* [HHb]';
    data_HHb=ft_selectdata(cfg, data_timelock_BL);
    % and rename labels such that they have the same name as HHb channels
    for i=1:length(data_HHb.label)
      tmp = strsplit(data_HHb.label{i});
      data_HHb.label{i}=tmp{1};
    end
    % then plot both on the lay-out
    cfg                   = [];
    cfg.channel = {'all', '-*a', '-*b', '-*c', '-*d'};
    cfg.showlabels        = 'yes';
    cfg.layout            = layout;
    cfg.showoutline       = 'yes';
    cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
    cfg.linecolor        = 'rb'; % O2Hb is showed in red, HHb in blue
    cfg.linewidth         = 2;
%     cfg.ylim              = [-0.2 0.4];
    cfg.ylim              = [-0.3 0.3];
    ft_multiplotER(cfg, data_O2Hb, data_HHb);title([type ' ' condition ' ' ID])
    
    % make topoplot
    cfg = [];
    cfg.channel = {'all', '-*a', '-*b', '-*c', '-*d'};
    data_O2Hb_long = ft_selectdata(cfg, data_O2Hb);
    cfg =[];
    cfg.channel = {'all', '-*a', '-*b', '-*c', '-*d'};
    cfg.parameter = 'avg';
    cfg.xlim = [0 3];
%     cfg.zlim = [-0.1 0.2];
    cfg.zlim = [-0.3 0.3];
    cfg.layout = layout;
    ft_topoplotER(cfg, data_O2Hb)
    title([type ' ' condition ' ' ID])
    end
end


