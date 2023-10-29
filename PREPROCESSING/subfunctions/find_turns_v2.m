function   [turns]=find_turns_v2(data, rot_tra, run, ID, varargin)
% FIND_TURNS_V2
% This function detects the 180 degree turns based on the pelvis
% position. The begin and end of the turn were then calculated as the 
% timepoints where the pelvis orientation crossed a regression line 
% with more than 5° (see also supplementary material).
% This function was partly based
% on:
% - Salarian A, et al. iTUG, a sensitive and reliable measure of mobility. (IEEE Trans Neural Syst Rehabil Eng) doi:10.1109/TNSRE.2010.2047606
% - Miller Koop M, et al. Quantifying turning behavior and gait in Parkinson's disease using mobile technology. (IBRO Rep.) doi:10.1016/j.ibror.2018.06.002
%
% Use as
%   [turns] = find_turns_v2(data, rot_tra, run, ID, varargin)
%
% INPUT:
%       data         = motion data in fieldtrip data structure (not
%       reframed yet)
%       rot_tra      = rotation and translation parameters outputted from
%       ROTATION_TRANSLATION
%       run          = run number of which the data was given
%       ID           = ID of which the data was given
% Additional options can be specified in key-value pairs and can be:
%       'vis'    true or false (default = true)
%
% OUTPUT
%       turns = table with the turn events (onset, duration, type,
%       value)
%
% DEPENDENCIES:
% - q2e.m (see: external): function to convert quaternions to euler angles
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample; % sampling frequency

if strcmp(ID, 'PD16') & run==2 % empty table (no turns)
  turns=table('Size', [0 4], 'VariableTypes', {'double', 'double', 'string', 'string'}, 'VariableNames', {'onset', 'duration', 'type', 'value'});
  return
end

%% select pelvis orientation data
if vis
    figure; hold on;
end

cfg=[];
cfg.channel={'seg_Pelvis_orientation_Q0','seg_Pelvis_orientation_Q1', 'seg_Pelvis_orientation_Q2', 'seg_Pelvis_orientation_Q3'};
ang_pelvis=ft_selectdata(cfg, data);

% convert orientations to euler angles
[x, y, z]=q2e(ang_pelvis.trial{1}(1,:), ang_pelvis.trial{1}(2,:),ang_pelvis.trial{1}(3,:),ang_pelvis.trial{1}(4,:)); % convert orientation to euler angles
orient=rad2deg([x;y;z]);

% recalibrate orientation
start_orient=nanmedian(orient(3,[1:fs*60]));
if strcmp(ID, 'PD17') & run==4
  start_orient = 70; % exception because of big drift at start of run
end
rot_tra.rot_angle;
turn_angle=orient(3,:)-start_orient;

% unwrap
turn_angle = unwrap(turn_angle, 340);

% low pass filter
turn_angle_lp=ft_preproc_lowpassfilter(turn_angle, fs, 5);
if strcmp(ID, 'PD57') & run==2
  % weird peak detected
  turn_angle_lp(17620:18115) = turn_angle_lp(17621);
elseif strcmp(ID, 'PD63') & run==4
  turn_angle_lp(11780:11880) = turn_angle_lp(11780);
end

if vis
  plot(turn_angle_lp)
end

%% find changepoints as the approximate begin and end of each turn
% find samples where orientation crosses 90 degrees, rounded at 10 degrees.
% Based on this, we know the number of turns
crossing=find(mod(round(turn_angle_lp, -1),180)==90);
crossing=[crossing(1) crossing(find([0 diff(crossing)]>600))];% if consecutive samples, use first sample
crossing=crossing(find(crossing>60*fs)); % e.g. PD90, run 5
nmb_chpt = length(crossing)*2;
if vis
  plot(crossing, turn_angle_lp(crossing), 'ro');
end

% find changepoints as the approximate begin and end of each turn
[icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', nmb_chpt, 'Statistic', 'linear');
% exceptions:
if strcmp(ID, 'PD57') % exception because very distorted turns with long FOG
  if run ==1
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 15, 'Statistic', 'linear');
  elseif run==2
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 17, 'Statistic', 'linear');
    icp = icp([1:3 5:10 12:13 15:17]);
  elseif run==4
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 13, 'Statistic', 'linear');
    icp = [icp length(turn_angle_lp)];
  end
elseif strcmp(ID, 'PD96')
  if run==1
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 17, 'Statistic', 'linear');
    icp = icp([1:3 5:12 14:15 17]);
  elseif run==5
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 17, 'Statistic', 'linear');
    icp = icp([1 3:12 14:15 17]);
  end
elseif strcmp(ID, 'PD17')
  if run==3
    [icp] = findchangepts(turn_angle_lp, 'MaxNumChanges', 15, 'Statistic', 'linear');
    turn_angle_lp(icp(10):icp(11))=2*turn_angle_lp(icp(10))-turn_angle_lp(icp(10):icp(11)); %weird signal change detected
    turn_angle_lp(icp(11)+1:end) = turn_angle_lp(icp(11)+1:end)-turn_angle_lp(icp(11)+1)+turn_angle_lp(icp(11));
    icp = icp([1:9 11:15]);
    if vis
    plot(turn_angle_lp)
    end
  end
end

% fit linear regression between the changepts
istart = [1 icp];
istop = [icp-1 length(turn_angle_lp)];
seg = struct();
for i = 1:length(icp)+1
  x = [istart(i):istop(i)];
  [P,S] = polyfit(x, turn_angle_lp(x), 1);
  [y_fit, delta] = polyval(P, x, S);
  if vis
  if i<=length(icp)
    xline(icp(i));
  end
  plot(x, y_fit, 'r-');
  plot(x, y_fit+2*delta, 'm--', x, y_fit-2*delta, 'm--');
  end
  % save for further processing
  seg(i).x = x;
  seg(i).y_fit = y_fit;
  seg(i).delta = delta;
end

%% define the exact begin & end of the turns
% these are define as the timepoint where the pelvis crossed the regression
% line with more than 5 degrees
turns=table('Size', [length(crossing) 4], 'VariableTypes', {'double', 'double', 'string', 'string'}, 'VariableNames', {'onset', 'duration', 'type', 'value'}); % create an empty table to store the information
idx=1; % number of the turn
margin = 5; % pelvis orientation crosses the regression line with more than 5 degrees.

for i=crossing
    % exception:
  if strcmp(ID, 'PD89') & run==1 & i==crossing(6)
    turns(idx,:)=[];
    % unexpected movement
    continue
  end
  % find the changepoint before and after the crossing
  chpt_before = find(icp<i,1, 'last');
  chpt_after = find(icp>=i, 1, 'first');
  if turn_angle_lp(icp(chpt_before))<turn_angle_lp(icp(chpt_after)) % turn to the left
    turns.value(idx)='turn_180_L';
    threshold = seg(chpt_before).y_fit(end)+margin;
    x_temp = [istart(chpt_before):istop(chpt_after)];
    inters_before = find(turn_angle_lp(x_temp)<threshold, 1, 'last');
    start_turn = x_temp(inters_before);
    threshold = seg(chpt_after+1).y_fit(1)-margin;
    x_temp = [istart(chpt_after):istop(chpt_after+1)];
    inters_after = find(turn_angle_lp(x_temp)>threshold, 1, 'first');
    stop_turn = x_temp(inters_after);
    % plot
    if vis
      plot([start_turn stop_turn], turn_angle_lp([start_turn stop_turn]), 'g*');
    end
  elseif turn_angle_lp(icp(chpt_before))>turn_angle_lp(icp(chpt_after)) % turn to the right
    turns.value(idx)='turn_180_R';
    threshold = seg(chpt_before).y_fit(end)-margin;
    x_temp = [istart(chpt_before):istop(chpt_after)];
    inters_before = find(turn_angle_lp(x_temp)>threshold, 1, 'last');
    start_turn = x_temp(inters_before);
    threshold = seg(chpt_after+1).y_fit(1)+margin;
    x_temp = [istart(chpt_after):istop(chpt_after+1)];
    inters_after = find(turn_angle_lp(x_temp)<threshold, 1, 'first');
    stop_turn = x_temp(inters_after);
    if isempty(stop_turn)
      loc_min=find(islocalmin(turn_angle_lp(seg(chpt_after+1).x), 'MinSeparation', fs/2) & turn_angle_lp(seg(chpt_after+1).x) < seg(chpt_after+1).y_fit + 3*seg(chpt_after+1).delta);
      stop_turn=seg(chpt_after+1).x(loc_min(1));% take first one
      warning('Using the other method')
    end
    if vis
      plot([start_turn stop_turn], turn_angle_lp([start_turn stop_turn]), 'g*');
    end
  end
  % check if turn is approx. 180 degrees
  if abs(turn_angle_lp(stop_turn)-turn_angle_lp(start_turn))<165 | abs(turn_angle_lp(stop_turn)-turn_angle_lp(start_turn))>190
    warning('Turn angle %d was %.0f degrees', idx, abs(turn_angle_lp(stop_turn)-turn_angle_lp(start_turn)))
  end
  % save in table
  turns.onset(idx)=start_turn/fs;
  turns.duration(idx)=(stop_turn-start_turn)/fs;
  turns.type(idx)='gait_event';
  idx=idx+1;
end

if vis
  title(sprintf('Right pelvis orientation around z-axis with detected turn events (sub-%s run-%02d)', ID, run))
  legend({'trunk orientation', 'pelvis orientation', 'crossing','change points', 'linear trend', 'delta 1', 'delta 2', 'start & stop of turn'});
end

