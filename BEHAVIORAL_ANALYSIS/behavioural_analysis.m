%% BEHAVIOURAL_ANALYSIS
% This script calculates and plots a variety of descriptional statistics
% based on the behavioral (motion) data. For instance, it can plot the walking
% speed during various conditions, calculates the average turn duration, or
% calculates at which point of the turn/doorway the FOG occured. 

%% Initialisation
root_dir='C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data';
figure_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\paper\figures\';
addpath(fullfile(root_dir,'scripts', 'final'))

%% Block averaging for motion data
clear data_timelock_temp data_timelock
variable = {'seg_Pelvis_velocity_X'}; % 'seg_Pelvis_acceleration_X'
group = {'PD', 'HC'};
type = {'normal'}; % FOG or normal
condition = {'turning walking'};
for t = length(type)
  for c = length(condition)
    for g = 1:length(group)
      switch group{g}
        case 'PD'
          ID={'PD06', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
            'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
            'PD76', 'PD77', 'PD88', 'PD90', 'PD96'};
        case 'HC'
          ID={'HC01', 'HC02', 'HC04', 'HC06', 'HC13', 'HC19', 'HC21', 'HC28', 'HC29', 'HC33', 'HC34', 'HC35'...
            'HC41', 'HC42', 'HC57', 'HC60', 'HC66', 'HC67', 'HC68', 'HC76', 'HC90', 'HC91'};
      end
      ID=sort(ID);
      for i=1:length(ID)
        load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_trialsmotion.mat', ID{i})));
        data_timelock_temp{g}{i} = timelock(trials_motion, ID{i}, 0, type{t}, condition{c}, false);
      end
      subjects=find(~cellfun(@isempty,data_timelock_temp{g}));% remove empty subjects
      data_timelock{g} = data_timelock_temp{g}(subjects);
      data_GA{g} = ft_timelockgrandaverage([], data_timelock{g}{:});
    end
    for v = 1:length(variable)
      figure; hold on;
      title([variable{v} ' ' condition{c}]);
      plotCI(data_GA{1}, variable{v}, length(data_GA{1}.cfg.previous), [248/255, 118/255, 109/255]);
      if ~strcmp(condition{c}, 'FOG')
        plotCI(data_GA{2},  variable{v}, length(data_GA{2}.cfg.previous), [0, 191/255, 196/255]);
      end
%       ylim([-0.5 1.2])
%       ylim([0 1.2]); ylabel('velocity (m/s)'); xlim([-10 10]); xlabel('time (s)')
      a = plot(1, nan, '-',  'Color', '#F8766D', 'linewidth', 2); b = plot(1,nan, '-', 'Color', '#00BFC4', 'linewidth', 2);
      legend([a, b], {'PD', 'HC'});
    end
  end
end

%% Statistics for motion data
% (requires first a block average of the Pelvis velocity/acceleration in
% the previous step)
% start & stop
% Is the average walking velocity different between the two groups?
cfg = [];
cfg.channel = 'seg_Pelvis_velocity_X';
cfg.latency = [7 10];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}+N{2}) = [ones(1,N{1}) 2*ones(1,N{2})];
stat = ft_timelockstatistics(cfg, data_timelock{1}{:}, data_timelock{2}{:})
% specific speed of each group:
cfg = [];
cfg.latency = [7 10];
cfg.channel = 'seg_Pelvis_velocity_X';
cfg.avgovertime = 'yes';
speed_PD = ft_selectdata(cfg, data_GA{1});
speed_HC = ft_selectdata(cfg, data_GA{2});
speed_PD.avg
speed_HC.avg
% Does one group starts quicker than the other?
cfg = [];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.latency = [0 3];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}+N{2}) = [ones(1,N{1}) 2*ones(1,N{2})];
stat = ft_timelockstatistics(cfg, data_timelock{1}{:}, data_timelock{2}{:})
% specific acceleration of each group:
cfg = [];
cfg.latency = [0 3];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.avgovertime = 'yes';
acc_PD = ft_selectdata(cfg, data_GA{1});
acc_HC = ft_selectdata(cfg, data_GA{2});
acc_PD.avg
acc_HC.avg
% Does one group stops quicker than the other?
cfg = [];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.latency = [-3 0];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}+N{2}) = [ones(1,N{1}) 2*ones(1,N{2})];
stat = ft_timelockstatistics(cfg, data_timelock{1}{:}, data_timelock{2}{:})
% specific acceleration of each group:
cfg = [];
cfg.latency = [-3 0];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.avgovertime = 'yes';
acc_PD = ft_selectdata(cfg, data_GA{1});
acc_HC = ft_selectdata(cfg, data_GA{2});
acc_PD.avg
acc_HC.avg

% turn
% Does one group slows more down before turning?
cfg = [];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.latency = [-2 0];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}+N{2}) = [ones(1,N{1}) 2*ones(1,N{2})];
stat = ft_timelockstatistics(cfg, data_timelock{1}{:}, data_timelock{2}{:})
% specific acceleration of each group:
cfg = [];
cfg.latency = [-2 0];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.avgovertime = 'yes';
acc_PD = ft_selectdata(cfg, data_GA{1});
acc_HC = ft_selectdata(cfg, data_GA{2});
acc_PD.avg
acc_HC.avg
% or better to calculate the min. acceleration? or where the decelerration
% starts? (similar for door)

% door
% Is the deceleration between the groups different?
cfg = [];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.latency = [-2 0];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}+N{2}) = [ones(1,N{1}) 2*ones(1,N{2})];
stat = ft_timelockstatistics(cfg, data_timelock{1}{:}, data_timelock{2}{:})
% specific acceleration of each group:
cfg = [];
cfg.latency = [-2 0];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.avgovertime = 'yes';
acc_PD = ft_selectdata(cfg, data_GA{1});
acc_HC = ft_selectdata(cfg, data_GA{2});
acc_PD.avg
acc_HC.avg
% Does each group slow down at the door?
cfg = [];
cfg.channel = 'seg_Pelvis_acceleration_X';
cfg.latency = [-2 0];
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'stats';
cfg.statistic = 'ttest';
N{1} = length(data_timelock{1});
N{2} = length(data_timelock{2});
cfg.design(1,1:N{1}) = [ones(1,N{1})];
stat_PD = ft_timelockstatistics(cfg, data_timelock{1}{:})
cfg.design = [];
cfg.design(1,1:N{2}) = [ones(1,N{2})];
stat_HC = ft_timelockstatistics(cfg, data_timelock{2}{:})

%% What was the duration of turns?
for g = 1:length(group)
  switch group{g}
    case 'PD'
      ID={'PD06', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
        'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
        'PD76', 'PD77', 'PD88', 'PD90', 'PD96'};
    case 'HC'
      ID={'HC01', 'HC02', 'HC04', 'HC06', 'HC13', 'HC19', 'HC21', 'HC28', 'HC29', 'HC33', 'HC34', 'HC35'...
        'HC41', 'HC42', 'HC57', 'HC60', 'HC66', 'HC67', 'HC68', 'HC76', 'HC90', 'HC91'};
  end
  ID=sort(ID);
  turn_duration{g} = nan(length(ID),1);
  for i=1:length(ID)
    load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_trialsmotion.mat', ID{i})));
    trials_turn = find(trials_motion.trialinfo(:,1) == 2 & (trials_motion.trialinfo(:,2) == 111|trials_motion.trialinfo(:,2) == 121));
    turn_duration{g}(i) = nanmean(trials_motion.trialinfo(trials_turn, 5));
  end
end
[H,P,CI,STATS] = ttest2(turn_duration{1}, turn_duration{2});
nanmean(turn_duration{1})
nanmean(turn_duration{2})

%% Where at the doorway did the patients freeze? Where at the turn?
ID={'PD06', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
  'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
  'PD76', 'PD77', 'PD88', 'PD90', 'PD96'};
door = [];
turn = [];
for s=1:length(ID)
  fprintf('====== sub-%s ======= \n', ID{s})
  load(fullfile(root_dir, 'processed', ['sub-' ID{s}], sprintf('sub-%s_trialsmotion.mat', ID{s})));
  sel_doorfog = (trials_motion.trialinfo(:,1)==1 & trials_motion.trialinfo(:,2)>20 & trials_motion.trialinfo(:,2)<25);
  sel_turnfog = (trials_motion.trialinfo(:,1)==1 & trials_motion.trialinfo(:,2)>100);
  door = [door; repmat(s, sum(sel_doorfog), 1) trials_motion.trialinfo(sel_doorfog, 4)]; % extract the distance trigger info
  sel_turnfog_during = sel_turnfog & trials_motion.trialinfo(:,4)>0; % FOG happening before the turn
  sel_turnfog_before = sel_turnfog & trials_motion.trialinfo(:,4)<=0; % FOG happening after the turn
  turn = [turn; repmat(s, sum(sel_turnfog_during), 1) trials_motion.trialinfo(sel_turnfog_during, 4); repmat(s, sum(sel_turnfog_before), 1) trials_motion.trialinfo(sel_turnfog_before, 4)]; 
end

% plot the door data
doorfreezers = unique(door(:,1));
n = length(doorfreezers);
door_median = nan(n, 1);
figure; hold on;
for s=1:length(doorfreezers)
  idx = find(door(:,1) == doorfreezers(s));
  plot(door(idx, 2)', n + rand(1,length(idx))/2, 'o');
  n = n-1;
  door_median(s) = median(door(idx, 2));
end
yticks([1:length(doorfreezers)]);
yticklabels(ID(flip(doorfreezers)));
xline(0, 'Label', 'door', 'LineWidth', 2);
xlabel('distance (m)')
set(gca, 'fontsize', 9)
set(gcf, 'Units', 'centimeters')
set(gcf,'Position',[1 1 9 9]);
filename = fullfile(figure_dir, 'FOGoccurence_plots', 'door');
% print(filename, '-depsc')
% saveas(gcf, [filename '.png'])

% calculate median & IQR
median(door_median)
quantile(door_median, [0.25 0.75])

% plot the turn data
turnfreezers = unique(turn(:,1));
n = length(turnfreezers);
turn_median = nan(n, 1);
figure; hold on;
for s=1:length(turnfreezers)
  idx = find(turn(:,1) == turnfreezers(s));
  plot(turn(idx, 2)', n + rand(1,length(idx))/2, 'o');
  n = n-1;
%   turn(idx(turn(idx,2)<0),2) = 0; % replace distance from turn by 0 (FOG at angle 0) (only for statistics)
  turn_median(s) = median(turn(idx, 2));
end
yticks([1:length(turnfreezers)])
yticklabels(ID(flip(turnfreezers)))
xline(0, 'Label', 'start turn', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
xline(1, 'Label', 'end turn',  'LineWidth', 2);
xlim([-0.2 1.1])
xlabel('turn angle (%)')
set(gca, 'fontsize', 9)
set(gcf, 'Units', 'centimeters')
set(gcf,'Position',[1 1 9 9]);
filename = fullfile(figure_dir, 'FOGoccurence_plots', 'turn');
% print(filename, '-depsc')
% saveas(gcf, [filename '.png'])

% calculate median & IQR
median(turn_median)*180
quantile(turn_median, [0.25 0.75])*180

%% Percentage of turns L/R that leaded to FOG
ID={'PD06', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
  'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
  'PD76', 'PD77', 'PD88', 'PD90', 'PD96'};
turns = zeros(length(ID), 2); % L - R
FOGs = zeros(length(ID), 2); % L - R
FOG_duration = zeros(length(ID), 2);
for s=1:length(ID)
  fprintf('===== %s =====\n', ID{s})
  load(fullfile(root_dir, 'processed', ['sub-' ID{s}], sprintf('sub-%s_run.mat', ID{s})));
  for r = 1:length(run)
    if ~isempty(run(r).data_nirs)
      run_turns = run(r).events(contains(run(r).events.value, '180'), :);
      turn_nrl = run_turns(strcmp(run_turns.type, 'gait_event'),:);
      turn_FOG = run_turns(strcmp(run_turns.type, 'FOG_Trigger'),:);
      turns(s,1) = turns(s,1) + sum(contains(turn_nrl.value, '_L_'));
      turns(s,2) = turns(s,2) + sum(contains(turn_nrl.value, '_R_'));
      FOGs(s,1) = FOGs(s,1) + sum(contains(turn_FOG.value, '_L_'));
      FOGs(s,2) = FOGs(s,2) + sum(contains(turn_FOG.value, '_R_'));
      FOG_duration(s,1) = FOG_duration(s,1) + sum(turn_FOG.duration(contains(turn_FOG.value, '_L_')));
      FOG_duration(s,2) = FOG_duration(s,2) + sum(turn_FOG.duration(contains(turn_FOG.value, '_R_')));
    end
  end
end
percentage = FOGs./turns*100
percentage_frozen = FOG_duration./turns
frozen_corrected = FOG_duration.*sum(turns,2)./turns
% save('LvRturns.mat', 'turns', 'FOGs', 'FOG_duration')
