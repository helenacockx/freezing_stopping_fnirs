%% FOG_EVENTS
% Analyzes the FOG events and provides descriptional statistics of the
% observed FOG events
%
% INPUT:
% - *run.mat files (contains all events)
% 
% OUTPUT:
% - descriptional statistics in the command window
% - plots of the descriptional statistics
% - FOG_events_all.mat: table of all events of all participants (can be
% reloaded for this function, so data does not need to be regenerated)

%% Initialization
root_dir='C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data';

ID={'PD06', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
  'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
  'PD76', 'PD77', 'PD88', 'PD90', 'PD96'};

ID=sort(ID);

%% load events
% load gait events of all participants
gait_events_all=[]; 
for i=1:length(ID)
  fprintf('----- Loading events of %s...  -----\n', ID{i})
  load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_run.mat', ID{i})));
  for r=1:length(run)
    if isempty(run(r).events)
      continue
    end
    run(r).events.ID= repmat(ID{i}, height(run(r).events), 1);
    gait_events_run=run(r).events(contains(run(r).events.type, 'gait_event'),:);
    gait_events_all=[gait_events_all; gait_events_run];
  end
end

% load FOG events and gait tasks of all participants
FOG_Trig_all=[]; FOG_Type_all=[]; gait_tasks_all=[];
for i=1:length(ID)
  fprintf('----- Loading events of %s...  -----\n', ID{i})
  load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_run.mat', ID{i})));
  for r=1:length(run)
    if isempty(run(r).events)
      continue
    end
    run(r).events.ID= repmat(ID{i}, height(run(r).events), 1);
    FOG_Trig_run=run(r).events(contains(run(r).events.type, 'FOG_Trigger'),:);
    FOG_Type_run=run(r).events(contains(run(r).events.type, 'FOG_Type'),:);
    gait_tasks=run(r).events(strcmp(run(r).events.type, 'gait_task'),:);
    FOG_Trig_all=[FOG_Trig_all; FOG_Trig_run];
    FOG_Type_all=[FOG_Type_all; FOG_Type_run];
    gait_tasks_all = [gait_tasks_all; gait_tasks];
  end
end
% save
FOG_events_all.FOG_Trig_all = FOG_Trig_all;
FOG_events_all.FOG_Type_all = FOG_Type_all;
FOG_events_all.gait_tasks_all = gait_tasks_all;
save(fullfile(root_dir, 'processed', 'FOG_events_all.mat'), 'FOG_events_all');

% (once loaded, we can also load this data from memory):
load(fullfile(root_dir, 'processed', 'FOG_events_all.mat'));
FOG_Trig_all = FOG_events_all.FOG_Trig_all;
FOG_Type_all = FOG_events_all.FOG_Type_all;
gait_tasks_all = FOG_events_all.gait_tasks_all;

%% generate general info about the FOG events
total_number_FOG=height(FOG_Trig_all)
total_number_FOG- sum(contains(FOG_Trig_all.value, 'standing')) % without the FOGs during standing
median_duration_FOG=median(FOG_Trig_all.duration) % in seconds
iqr_duration_FOG=prctile(FOG_Trig_all.duration, [25 75])
min_duration_FOG=min(FOG_Trig_all.duration)
max_duration_FOG=max(FOG_Trig_all.duration)

% generate information about the FOG events per participant
table_FOGmeanduration = grpstats(FOG_Trig_all, "ID", ["mean", "median", "min", "max"], "DataVars", "duration", "VarNames", ["ID", "FOGcount", "mean_FOGduration", "median_FOGduration","min_FOGduration", "max_FOGduration"]);
table_FOGtotduration = grpstats(FOG_Trig_all, "ID", @(x)sum(x,1), "DataVars", "duration", "VarNames", ["ID", "FOGcount", "total_FOGduration"]);
table_gaittask = grpstats(gait_tasks_all, "ID", @(x)sum(x,1), "DataVars", "duration", "VarNames", ["ID", "task_count", "total_taskduration"]);
temp_t = outerjoin(table_FOGmeanduration, table_FOGtotduration,'Keys', {'ID', 'FOGcount'}, 'MergeKeys', true);
part_table = outerjoin(temp_t, table_gaittask,'Keys', 'ID', 'MergeKeys', true);
[m,n] = find(isnan(part_table{:,[2:end]}));
part_table{m,n+1} = 0;
part_table.TF = part_table.total_FOGduration./part_table.total_taskduration*100 % TF = percentage time frozen

% summary of the information about the FOG events per participant
num_noFOG = sum(part_table.FOGcount == 0)
median_number_FOG_part = median(part_table.FOGcount)
iqr_number_FOG_part = prctile(part_table.FOGcount, [25 75])
median_duration_FOG_part=median(part_table.median_FOGduration) % in seconds
iqr_duration_FOG_part=prctile(part_table.median_FOGduration, [25 75])
mean_TF = mean(part_table.TF)
median_TF = median(part_table.TF)
iqr_TF = prctile(part_table.TF, [25 75])


%% boxplots of FOG durations
figure; boxplot(FOG_Trig_all.duration, FOG_Trig_all.ID, 'PlotStyle', 'traditional', 'Colors','k', 'Symbol', 'k.', 'DataLim',[0 75], 'width', 0.8)
 xlabel('participant'); ylabel('FOG duration (s)')
 title('FOG duration by participant')
%  set(gca, 'Fontsize', 40);
%  set(findobj(gca,'type','line'),'linew',3)
% saveas(gcf, fullfile(fig_dir, 'FOG_variables','FOGduration_byPatient.eps'))

%% make bar charts showing differences between each patient
% by trigger
n=length(ID);
nmb_trig=nan(n,5); dur_trig=nan(n,5); nmb_trig_LR=nan(n,2); dur_trig_LR=nan(n,2);
for i=1:length(ID)
  total_dur=sum(gait_tasks_all.duration(strcmp(gait_tasks_all.ID, ID(i))));
  idx_turn= find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, '180'));
  idx_turn_left = find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, '180') & contains(FOG_Trig_all.value, '_L_'));
  idx_turn_right = find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, '180') & contains(FOG_Trig_all.value, '_R_'));
  idx_door=find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, 'Doorway')); 
  idx_tar=find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, 'Target'));
  idx_sh=find(strcmp(FOG_Trig_all.ID, ID(i)) & contains(FOG_Trig_all.value, 'SH'));
  idx_other=find(strcmp(FOG_Trig_all.ID, ID(i)) & strcmp(FOG_Trig_all.value, 'FOG'));
  nmb_trig(i,1)=length(idx_turn);
  nmb_trig(i,2)=length(idx_door);
  nmb_trig(i,3)=length(idx_tar);
  nmb_trig(i,4)=length(idx_sh);
  nmb_trig(i,5)=length(idx_other);
  nmb_trig_LR(i,1)=length(idx_turn_left);
  nmb_trig_LR(i,2)=length(idx_turn_right);
  dur_trig(i,1)=sum(FOG_Trig_all.duration(idx_turn))/total_dur*100;
  dur_trig(i,2)=sum(FOG_Trig_all.duration(idx_door))/total_dur*100;
  dur_trig(i,3)=sum(FOG_Trig_all.duration(idx_tar))/total_dur*100;
  dur_trig(i,4)=sum(FOG_Trig_all.duration(idx_sh))/total_dur*100;
  dur_trig(i,5)=sum(FOG_Trig_all.duration(idx_other))/total_dur*100;
  dur_trig_LR(i,1)=sum(FOG_Trig_all.duration(idx_turn_left))/total_dur*100;
  dur_trig_LR(i,2)=sum(FOG_Trig_all.duration(idx_turn_right))/total_dur*100;
end
% for number of freezes
figure;
bar(nmb_trig, 'stacked'); % or: bar(nmb_trig);
ylim([0 100]); xticks([1:length(ID)]); xticklabels(ID);
xlabel('patient'); ylabel('number of freezes')
legend({'FOG turn', 'FOG narrow passage','FOG target', 'FOG starting hesitation', 'FOG others'})
title('Number of FOGs by patient and by trigger')
% saveas(gcf, fullfile(fig_dir, 'FOG_variables','NumberFOGs_byPatient_byTrigger.jpg'))
% for total time frozen
figure;
bar(dur_trig, 'stacked'); % or: bar(nmb_trig, 'stacked');
ylim([0 50]); xticks([1:length(ID)]); xticklabels(ID);
xlabel('patient'); ylabel('percentage of total time frozen')
legend({'FOG turn', 'FOG narrow passage','FOG target','FOG starting hesitation', 'FOG others'})
title('Percentage of total time frozen by patient and by trigger')
% % saveas(gcf, fullfile(fig_dir, 'FOG_variables','PercentageFOGs_byPatient_byTrigger.jpg'))
% in numbers
prct = nmb_trig./sum(nmb_trig,2)*100;
prct(isnan(prct))=0;
mean(prct)
sum(nmb_trig)./sum(sum(nmb_trig))*100 % percentage of all the FOG episodes
sum((nmb_trig)~=0) % number of participants with this type of freezing

% by type
nmb_type=nan(n,2); dur_type=nan(n,2);
for i=1:n
  total_dur=sum(gait_tasks_all.duration(strcmp(gait_tasks_all.ID, ID(i))));
  idx_tremb= find(strcmp(FOG_Type_all.ID, ID(i)) & contains(FOG_Type_all.value, 'Trembling'));
  idx_akin=find(strcmp(FOG_Type_all.ID, ID(i)) & contains(FOG_Type_all.value, 'Akinesia'));
  idx_shuf=find(strcmp(FOG_Type_all.ID, ID(i)) & contains(FOG_Type_all.value, 'Shuffling'));
  nmb_type(i,1)=length(idx_tremb) + length(idx_shuf);
  nmb_type(i,2)=length(idx_akin);
  dur_type(i,1)=sum(FOG_Type_all.duration(idx_tremb))/total_dur*100 + sum(FOG_Type_all.duration(idx_shuf))/total_dur*100;
  dur_type(i,2)=sum(FOG_Type_all.duration(idx_akin))/total_dur*100;
end
% for number of freezes
figure;
bar(nmb_type, 'stacked'); % or: bar(nmb_trig);
ylim([0 100]); xticks([1:n]);xticklabels(ID);
xlabel('patient'); ylabel('number of freezes')
legend({'trembling/shuffling', 'akinesia'})
title('Number of FOGs by patient and by type')
% saveas(gcf, fullfile(fig_dir, 'FOG_variables','NumberFOGs_byPatient_byType.jpg'))
% for total time frozen
figure;
bar(dur_type, 'stacked'); % or: bar(nmb_trig, 'stacked');FOG
ylim([0 50]); xticks([1:n]);xticklabels(ID);
xlabel('patient'); ylabel('percentage of total time frozen')
legend({'trembling/shuffling',  'akinesia'})
title('Percentage of total time frozen by patient and by type')
% saveas(gcf, fullfile(fig_dir, 'FOG_variables','PercentageFOGs_byPatient_byType.jpg'))
% in numbers
prct = nmb_type./sum(nmb_type,2)*100;
prct(isnan(prct))=0;
mean(prct)
sum(nmb_type)./sum(sum(nmb_type))*100 % percentage of all the FOG episodes
sum((nmb_type)~=0) % number of participants with this type of freezing

% make also pie charts
figure; subplot(1,2,1)
pie(sum(nmb_trig), {'FOG turn', 'FOG door', 'FOG target', 'FOG starting hesitation', 'FOG others'}); 
subplot(1,2,2)
pie(sum(nmb_type), {'trembling/shuffling', 'akinesia'}); 
% split in left and right turn
pie_data = [sum(nmb_trig_LR) sum(nmb_trig(:,2:end))];
figure;
pie(pie_data, {'FOG turn left' 'FOG turn right', 'FOG door', 'FOG target', 'FOG starting hesitation', 'FOG others'}); 
% pie for FOG severity
figure;
pie(sum(dur_trig), {'FOG turn', 'FOG door', 'FOG target', 'FOG starting hesitation', 'FOG others'}); 
pie_data = [sum(dur_trig_LR) sum(dur_trig(:,2:end))];
figure;
pie(pie_data, {'FOG turn left' 'FOG turn right', 'FOG door', 'FOG target', 'FOG starting hesitation', 'FOG others'}); 
