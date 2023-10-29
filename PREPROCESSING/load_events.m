function [run] = load_events(run, events, ID)
%% LOAD_EVENTS
% Loads the (ELAN) events for the given participant and stores it in the run
% structure. It also redefines some of the motion events based on the
% motion data by using a better function (v2) than was used for the import events for
% ELAN.
%
% INPUT:
% - run: matlab structure containing the raw nirs & motion data of each run
% - ID: participant-ID of the run that is currently loaded
% - events: table with the sync events for this participant of the whole
% session, with t = 0 as the start of the first run.
% - *annotations-Helena.txt: exported ELAN annotations performed by H.C. for this
% participant.
% - *annotations-Yuli.txt: exported ELAN annotations performed by Y.A.F.R.
% for this participant.
%
% OUTPUT:
% - run: same as input, but field 'events' added for each run, containing
% all motion and FOG events, synchronized to the start of the run (t = 0).
% 
% DEPENDENCIES: 
% - rotation_translation.m:
% - find_turns_v2.m
% - find_startstop_v1.m
% - find_startstop_v2.m
%
%%
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

% directories containing the ELAN annotations
annotations_H_dir = '\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\Helena';
annotations_agreed_dir = '\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data\processed\annotations\agreed';

%% load gait_event, gait_task, and unexpected movements from annotations of Helena
fprintf('----- Loading events of %s...  -----\n', ID)
filename = fullfile(annotations_H_dir, ['sub-' ID '_annotations-Helena.txt']);
opts = detectImportOptions(filename);
annotations_H = readtable(filename, opts);
annotations_H_long = stack(annotations_H, 3:size(annotations_H,2), 'NewDataVariableName', 'value', 'IndexVariableName', 'type' );
annotations_H_long = annotations_H_long(~ismissing(annotations_H_long.value) & ~ismissing(annotations_H_long.value, ''),:); % remove empty rows
annotations_H_long.type = cellstr(annotations_H_long.type); % convert categorical to string array
annotations_H_long.value = cellstr(annotations_H_long.value); % convert charactero to string array
annotations_H_long.EndTime_Ss_msec = annotations_H_long.EndTime_Ss_msec - annotations_H_long.BeginTime_Ss_msec;
annotations_H_long.Properties.VariableNames = {'onset', 'duration', 'type', 'value'};

% select relevant events
events_H = annotations_H_long(contains(annotations_H_long.type(:), {'gait_task', 'Unexpected_movements', 'gait_event'}) & ~strcmp(annotations_H_long.value(:), 'first_heel'),:);

% make teh following events 0 s duration: start_walking, stop_walking, doorway
zero_events = find(contains(events_H.value, {'start_walking', 'stop_walking', 'doorway'}));
events_H.duration(zero_events) = 0;

%% load agreed FOG_events
try
  filename = fullfile(annotations_agreed_dir, ['sub-' ID '_annotations-agreed.txt']);
  opts = detectImportOptions(filename);
  annotations_agr =  readtable(filename, opts);
  annotations_agr = annotations_agr(:, {'BeginTime_Ss_msec', 'EndTime_Ss_msec', 'FOG_agreed_Trigger',  'FOG_agreed_Type', 'FOG_disagreed_Trigger', 'FOG_disagreed_Type'});
  annotations_agr = annotations_agr(~ismissing(annotations_agr.FOG_agreed_Trigger, '') | ~ismissing(annotations_agr.FOG_disagreed_Type, ''),:); % remove empty rows

  % combine FOG_annotations that are neighbouring
  idx = find(annotations_agr.BeginTime_Ss_msec(2:end)-annotations_agr.EndTime_Ss_msec(1:end-1)<0.1);
  for i=idx'
    trigger = unique([annotations_agr.FOG_agreed_Trigger(i); annotations_agr.FOG_disagreed_Trigger(i); annotations_agr.FOG_agreed_Trigger(i+1); annotations_agr.FOG_disagreed_Trigger(i+1)]);
    trigger = trigger(~ismissing(trigger, ''));
    type = unique([annotations_agr.FOG_agreed_Type(i); annotations_agr.FOG_disagreed_Type(i); annotations_agr.FOG_agreed_Type(i+1); annotations_agr.FOG_disagreed_Type(i+1)]);
    type = type(~ismissing(type, ''));
    if length(trigger)>1 | length(type)>1
      % keep the FOG events this way
      idx(idx==i) = [];
      continue
    end
    % combine two annotations
    annotations_agr.EndTime_Ss_msec(i) = annotations_agr.EndTime_Ss_msec(i+1);
  end
  % delete spurious annotations
  annotations_agr(idx+1,:) = [];

  % change format
  annotations_agr_long = stack(annotations_agr, 3:size(annotations_agr,2), 'NewDataVariableName', 'value', 'IndexVariableName', 'type');
  annotations_agr_long = annotations_agr_long(~ismissing(annotations_agr_long.value(:), ''),:); % remove empty rows
  annotations_agr_long.type = cellstr(annotations_agr_long.type); % convert categorical to string array
  annotations_agr_long.value = cellstr(annotations_agr_long.value);
  annotations_agr_long.EndTime_Ss_msec = annotations_agr_long.EndTime_Ss_msec - annotations_agr_long.BeginTime_Ss_msec;
  annotations_agr_long.Properties.VariableNames = {'onset', 'duration', 'type', 'value'};

  % rename type
  events_FOG = annotations_agr_long;
  events_FOG.type(contains(events_FOG.type, 'Trigger')) = {'FOG_Trigger'};
  events_FOG.type(contains(events_FOG.type, 'Type')) = {'FOG_Type'};

catch
  % use annotations_H (PD17: no FOG)
  events_FOG = annotations_H_long(contains(annotations_H_long.type, 'FOG'),:);
end


%% combine events
events_all = [events_H; events_FOG];
events_all = sortrows(events_all);
if strcmp(ID, 'PD76')
  warning('acq-mobile was started 25.361 seconds too late. Annotations were shifted to sync them in ELAN. We will shift them back now.')
  events_all.onset(2:end)=[events_all.onset(2:end)] + 25.361; % not the first gait_task
  events_all.duration(1) = events_all.duration(1) + 25.361;
elseif strcmp(ID, 'HC28')
  warning('acq-mobile was started 80.639 seconds too late. Annotations were shifted to sync them in ELAN. We will shift them back now.')
  events_all.onset(3:end)=[events_all.onset(3:end)] + 80.639; % not the first gait_task & FOG course
  events_all.onset(1) = 62.822; % recover first gait event
  events_all.duration(2) = events_all.duration(2) + 80.639; % add delay to FOG_course
end


%% split events in runs and redefine some of the motion events (based on the motion data)
% remove the lsl FOG events (button press during experiment when FOG was
% observed)
events = events(~strcmp(events.value(:), 'FOG'),:);
events_all = sortrows([events; events_all]);

% split in runs
start_runs=events.onset(find(strcmp(events.value, 'start_run')));
stop_runs=events.onset(find(strcmp(events.value, 'stop_run')));
if strcmp(ID, 'PD06') % exception: events start from run-02 on (run-01 is empty)
  start_runs = [nan; start_runs];
  stop_runs = [nan; stop_runs];
elseif strcmp(ID, 'HC06') % run 1-4 are empty
  start_runs = [nan(4,1); start_runs];
  stop_runs = [nan(4,1); stop_runs];
end
numb_run = length(run);
for r=1:numb_run % run over runs
  if isempty(run(r).events) | isempty(run(r).data_motion)
    continue
  end
  % select events that have onsets falling within the time window of the
  % current run (with 3 s margin) and sync to the start of the run
  run_events=events_all((events_all.onset(:)>=start_runs(r)-3 & events_all.onset(:)<=stop_runs(r)+3),:);
  run_events.onset=run_events.onset(:)-start_runs(r); % t = 0 is onset of the run
  
  % trim data & events based on valid window (see run.info)
  rmvd_startblocks = [];
  rmvd_stopblocks = [];
  if strcmp(run(r).info.complete, 'no') 
    if run(r).info.validwindow(1)==0
      cfg=[];
      cfg.latency=run(r).info.validwindow;
      run(r).data_motion.data_raw=ft_selectdata(cfg, run(r).data_motion.data_raw);
      gait_task=run_events(strcmp(run_events.type, 'gait_task'),:); % adapt the gait task according to the valid window
      gait_task.onset = run(r).info.validwindow(1);
      gait_task.duration = run(r).info.validwindow(2)-gait_task.onset;
      idx_keep=find((run_events.onset>=cfg.latency(1) & run_events.onset<cfg.latency(2) & ~strcmp(run_events.type, 'gait_task')));
      rmvd_startblocks = run_events(strcmp(run_events.value, 'start_block') & (run_events.onset <= cfg.latency(1) | run_events.onset>=cfg.latency(2)),:);
      rmvd_stopblocks = run_events(strcmp(run_events.value, 'stop_block') & (run_events.onset <= cfg.latency(1) | run_events.onset>=cfg.latency(2)),:);
      run_events=run_events(idx_keep, :);
      run_events=[gait_task; run_events];
    end
  end
  
  %% redefine turns based on the new version of the function
  [rot_tra] = rotation_translation(run(r).data_motion.data_raw, r, ID, 'visualize', false);
  [turns]=find_turns_v2(run(r).data_motion.data_raw, rot_tra, r, ID, 'visualize', false);
  idx = find(contains(run_events.value, 'turn') & contains(run_events.type, 'gait_event'));
  if all(strcmp(run_events.value(idx), turns.value))
    run_events.onset(idx) = turns.onset(:);
    run_events.duration(idx) = turns.duration(:);
  else
    error
  end
  
  %% redefine starts and stops based on the new cutoff of 0.1 m/s
  % first check whether the start/stops were shifted during the annotation
  idx = find(contains(run_events.value, {'start', 'stop'}) & contains(run_events.type, 'gait_event'));
  [startstop_temp]=find_startstop_v1(run(r).data_motion.data_raw, run_events, rot_tra, r, ID, 'visualize', false);
  startstop_temp = sortrows(startstop_temp, 'onset');
  if strcmp(ID, 'PD16') & r==2
    % remove the last stop event
    startstop_temp = startstop_temp(1:end-1,:);
  elseif strcmp(ID, 'PD57') & r==4
    startstop_temp(end+1,:) = run_events(idx(end),:);
  end
  if any(abs(startstop_temp.onset - run_events.onset(idx)) > 0.1)
    shifted = find(abs(startstop_temp.onset - run_events.onset(idx)) > 0.1);
  else
    shifted = [];
  end
  % calculate the new start & stops
  [startstop]=find_startstop_v2(run(r).data_motion.data_raw, run_events, r, ID, 'visualize', false);
  startstop = sortrows(startstop, 'onset');
  if strcmp(ID, 'PD16') & r==2
    % remove the last stop event
    startstop = startstop(1:end-1,:);
  elseif strcmp(ID, 'PD57') & r==4
    startstop(end+1,:) = run_events(idx(end),:);
  end
  if ~isempty(shifted)
    % replace the start/stop event by the shifted annotation
    startstop(shifted,:) = run_events(idx(shifted),:);
  end
  if any(abs(startstop.onset - run_events.onset(idx)) > 1.5)
    warning('Discrepancy between the old and new event.')
    startstop.onset - run_events.onset(idx)
    run_events.onset(idx) = startstop.onset(:);
  elseif all(strcmp(run_events.value(idx), startstop.value))
    run_events.onset(idx) = startstop.onset(:);
  else
    error
  end
  
  %% add extra info to events with suffixes
  % add suffix for walking/standing condition
  start_blocks = run_events(strcmp(run_events.value, 'start_block'),:);
  start_blocks = sortrows([start_blocks; rmvd_startblocks]); % add removed start_blocks
  stop_blocks = run_events(strcmp(run_events.value, 'stop_block'),:);
  stop_blocks = sortrows([stop_blocks; rmvd_stopblocks]); % add_removed stop_blocks
  if height(start_blocks) ~= height(stop_blocks)
    error
  end
  idx2change = find(contains(run_events.value, {'180', 'doorway', '_SH'}, 'IgnoreCase', true));
  for i=idx2change'
    if any(run_events.onset(i)>[start_blocks.onset] & run_events.onset(i)<[stop_blocks.onset]) % during th blocks, so during walking
      run_events.value{i} = [run_events.value{i} '_walking'];
    else % in between the blocks, so during standing
      run_events.value{i} = [run_events.value{i} '_standing'];
    end
  end
  
  % add suffix to stop event to show whether the stop was in front of the
  % door/square (turn)
  idx2change = find(strcmp(run_events.value, 'stop_walking'));
  idxgaitevents = find(strcmp(run_events.type, 'gait_event'));
  for i=idx2change'
    idxprecedent = idxgaitevents(find(idxgaitevents<i, 1, 'last'));% find preceding gait_event
    if contains(run_events.value(idxprecedent), 'doorway')
      run_events.value{i} = [run_events.value{i} '_turn']; % if the precedent was a door, the stop is at a turn and vv
    elseif contains(run_events.value(idxprecedent), 'turn')
      run_events.value{i} = [run_events.value{i} '_door'];
    else
      error
    end
  end  
  
  % save events
  run(r).events = sortrows(run_events);
end
