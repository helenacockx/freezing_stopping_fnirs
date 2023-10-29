function   [startstop]=find_startstop_v2(data, events, run, ID, varargin)
% FIND_STARTSTOP_V2
% This function detects the timepoints when a participant starts and stops
% when walking. This was defined as an increase/decrease in horizontal
% pelvis velocity over a threshold of 0.1 m/s after a start/stop signal was given.
%
% Use as
%   [startstop]=find_startstop_v2(data, events, run, ID, varargin)
%
% INPUT:
%       data         = motion data in fieldtrip data structure (not
%       reframed yet)
%       events       = run events containing the start_block and stop_block
%       sync events
%       rot_tra      = rotation and translation parameters outputted from
%       ROTATION_TRANSLATION
%       run          = run number of which the data was given
%       ID           = ID of which the data was given
% Additional options can be specified in key-value pairs and can be:
%       'vis'    true or false (default = true)
%
% OUTPUT
%       startstop = table with the start and stop events (onset, duration, type,
%       value)
%
% DEPENDENCIES:/
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample;
%% select pelvis segment & transfer from global to local reference frame
cfg = [];
cfg.channel = find(contains(data.label, 'seg_Pelvis_velocity'));
data_pelvis = ft_selectdata(cfg, data);
cfg.channel = find(contains(data.label, 'seg_Pelvis_orientation'));
data_orient = ft_selectdata(cfg, data);
rotm = quat2rotm(data_orient.trial{1}');
pelvis_new = nan(3, length(data_pelvis.time{1}));
for i = 1:length(data_pelvis.time{1})
  pelvis_new(:,i) = rotm(:,:,i)'*data_pelvis.trial{1}(:,i);
end
data_pelvis.trial{1} = pelvis_new;


% plot
if vis
  figure; a=plot(data_pelvis.trial{1}(1,:)); hold on;
end

% low pass
cfg=[];
cfg.lpfilter='yes';
cfg.lpfreq=0.5;
% cfg.lpfreq=6;
data_pelvis_lp=ft_preprocessing(cfg, data_pelvis);

% combine X & Y for horizontal velocity
data_pelvis_lp.trial{1}(4,:) = sqrt(data_pelvis_lp.trial{1}(1,:).^2+data_pelvis_lp.trial{1}(2,:).^2);
data_pelvis_lp.label{4} = 'seg_Pelvis_velocity_horizontal';
if vis
  b=plot(data_pelvis_lp.trial{1}(4,:));
end

%% find start and stop walking events
% find start_walking as first point after start_block where velocity
% > threshold of 0.1 m/s
threshold = 0.1;
start_block=events(strcmp(events.value, 'start_block'),:);
onset=[]; value={}; idx=1;
for k=1:height(start_block)
  start_walking=find(abs(data_pelvis_lp.trial{1}(4,round(start_block.onset(k)*fs):end))>threshold,1,'first');
  onset(idx)=start_block.onset(k)+start_walking/fs;
  value{idx}='start_walking';
  if vis
    c=plot(onset(idx)*fs, data_pelvis_lp.trial{1}(4,round(onset(idx)*fs)), 'go');
  end
  idx=idx+1;
end

% find stop_walking as last point before stop_block where velocity <0.1
stop_block=events(strcmp(events.value, 'stop_block'),:);
for k=1:height(stop_block)
  if round(stop_block.onset(k)*fs+fs)>length(data_pelvis_lp.time{1})
    continue
  end
  stop_walking=find(abs(data_pelvis_lp.trial{1}(4,1:round(stop_block.onset(k)*fs+fs)))>threshold,1, 'last'); % stop_block button was sometime pressed earlier
  onset(idx)=stop_walking/fs;
  value{idx}='stop_walking';
  if vis
    d=plot(onset(idx)*fs, data_pelvis_lp.trial{1}(4,round(onset(idx)*fs)), 'ro');
  end
  idx=idx+1;
end

% create table
onset=onset';
value=value';
duration=zeros(size(onset));
type=repmat({'gait_event'}, size(onset));
startstop=table(onset, duration, type, value);

% vis
if vis
  title(sprintf('Pelvis velocity along the x-axis with detected start and stop events (sub-%s run-%02d)', ID, run))
  try
    legend([a,b,c,d], {'original', 'filtered', 'start walking', 'stop walking'});
  end
end
