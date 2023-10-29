function [data_epoch_nirs, data_epoch_motion] = define_trials(run, ID, vis, z_score)
% DEFINE_TRIALS
% This function segments the data into trials and adds info about each trial
%
% INPUT:
% - run: matlab structure containing the raw nirs & motion data of each run
% - ID: participant-ID of the run that is currently loaded
% - vis: visualize the trials? (1/0)
% - z_score: use z-scored data or not? (1/0)
%
% OUTPUT:
% - data_epoch_nirs: the fNIRS data segmented into trials (fieldtrip
% structure)
% - data_epoch_motion: the motion data segmented into trials
% (fieldtrip structure_
% The trial structures also contains trialinfo containing the following
% elements:
% - column 1 = type (FOG/nrl_gait)
%     1 = FOG
%     2 = normal gait event
% - column 2 = value (e.g. 180_R_walking/180_R_standing, doorway_walking,
% doorway_standing, ...)
%     1 = turn
%       11 = turn_R
%         111 = turn_R_walking
%         112 = turn_R_standing
%       12 = turn_L
%         121 = turn_L_walking
%         122 = turn_L_standing
%     2 = doorway
%       21 = doorway_walking
%       22 = doorway_standing
%     3 = FOG_Target/stop_walking
%       3 = FOG_target
%       3.1 = stop_walking_turn
%       3.2 = stop_walking_door
%     4 = FOG_SH/start_walking
%       4 = FOG_SH
%       4.1 = start_walking
%       4.2 = start_signal (start_block)
%     5 = FOG_others
% - column 3 = offset_trigger (if FOG, how long before/after the trigger
% occured the FOG? in sec., negative = FOG occured before trigger)
% - column 4 = distance_trigger (if FOG, how far before/after the trigger
% occured the FOG? if FOG turn, at what turn angle occured the FOG?)
% - column 5 = trial_duration
% - column 6 = run
%
% DEPENDENCIES:
% - quat2rotm.m
% - table2FT.m

run=load_run_info(run, ID);
tmp_nirs=struct([]);
tmp_motion=struct([]);

if strcmp(ID, 'PD77') | strcmp(ID, 'PD16') | strcmp(ID, 'PD90')
    margin_trig = 15; % turn was more than 10 sec away from FOG_target
else
    margin_trig = 10;
end
margin_FOG = 10;
excl_FOG = 0; excl_gait = 0;
offset_trigger_all=[];
distance_trigger_all =[];

for r=1:length(run)
    if isempty(run(r).data_nirs) | isempty(run(r).data_motion)
        continue
    end
    FOG_events = run(r).events(strcmp(run(r).events.type, 'FOG_Trigger'),:);
    gait_events = run(r).events(strcmp(run(r).events.type, 'gait_event'),:);
    sf = run(r).data_nirs.data_raw.fsample;
    offset = 10*sf;
    
    %% Define FOG events
    % exclude FOG events:
    % ...remove FOG_events (with 5-sec margin) that do not fall within the valid window
    idx = find(~([FOG_events.onset]-5>run(r).info.validwindow(1) & [FOG_events.onset+FOG_events.duration]+5 < run(r).info.validwindow(2)));
    FOG_events = FOG_events(setdiff(1:end, idx),:);
    excl_FOG = excl_FOG + length(idx);
    % ...remove FOG_events that fall within unexpected movements (5-sec
    % margin)
    if any(strcmp(run(r).events.type, 'Unexpected_movements'))
        unexp_mov = run(r).events(strcmp(run(r).events.type, 'Unexpected_movements'),:);
        for i=1:height(unexp_mov)
            idx = overlappingevt(FOG_events, unexp_mov.onset(i)-5, unexp_mov.onset(i)+unexp_mov.duration(i)+5);
            FOG_events = FOG_events(setdiff(1:end, idx),:);
            excl_FOG = excl_FOG + length(idx);
        end
    end
    % ...exclude FOG events that are preceded by another FOG event within 10
    % sec.
    idx = find([FOG_events.onset(2:end)-FOG_events.onset(1:end-1)-FOG_events.duration(1:end-1)] < margin_FOG);
    if ~isempty(idx)
        FOG_events = FOG_events(setdiff(1:end, idx+1),:);
        excl_FOG = excl_FOG + length(idx);
    end
    
    % define trials of the FOG events
    FOG_trl = [];
    FOG_eventsFT = table2FT(FOG_events, sf);
    trlbegin = [FOG_eventsFT(:).sample] - offset;
    trlend = [FOG_eventsFT(:).sample] + offset;
    type = ones(size(trlbegin));
    % add value info + offset_trigger
    value = nan(size(trlbegin)); offset_trigger = nan(size(trlbegin)); distance_trigger = nan(size(trlbegin)); trial_duration = nan(size(trlbegin));
    trigger_event_sample = nan(size(trlbegin)); % also store the triggering event
    for i=1:length(FOG_eventsFT)
        % find close gait_event
        idx = overlappingevt(gait_events, FOG_events.onset(i)-margin_trig, FOG_events.onset(i)+FOG_events.duration(i)+margin_trig);
        trigger_event = gait_events(idx,:);
        if isempty(trigger_event)
            if strcmp(ID, 'PD16') & r==2
                trigger_event = gait_events(2,:);
            else
                error('No close trigger event found for this FOG event')
            end
        end
        [value(i), trigger_event] = value_lookup(FOG_eventsFT(i).value, trigger_event);
        if height(trigger_event)>1
            % if multiple close events are detected, choose the closest one
            [M, idx]=min(abs(FOG_events.onset(i)-[trigger_event.onset]));
            trigger_event=trigger_event(idx,:);
        end
        if value(i) ==5
            offset_trigger(i) = nan;
        else
            offset_trigger(i) = FOG_events.onset(i) - trigger_event.onset;
        end
        
        trigger_eventFT = table2FT(trigger_event, sf);
        % for turn FOG, define turn angle of occurence of FOG
        if value(i)>100
            angle = calculate_turnangle(run(r).data_motion.data_raw, r, ID, trigger_eventFT, FOG_eventsFT(i));
            distance_trigger(i) = angle;
            % for doorway FOG, define the distance from the door at
            % the occurence of the FOG
        elseif value(i)>20
            distance_door = calculate_doordistance(run(r).data_motion.data_raw, r, ID, trigger_eventFT, FOG_eventsFT(i));
            distance_trigger(i) = distance_door;
        end
        
        % also provide information about the duration of the trial
        trial_duration(i) = FOG_events.duration(i);
        
    end
    FOG_trl    = [trlbegin' trlend' -1*offset*ones(size(trlbegin))' type' value' offset_trigger' distance_trigger' trial_duration' repmat(r, size(trlbegin'))];
    offset_trigger_all = [offset_trigger_all; value' offset_trigger'];
    distance_trigger_all = [distance_trigger_all; value' distance_trigger'];
    
    
    %% define normal gait events
    % = turn/door/start/stop without FOG in a 10-sec. margin
    % add start_block events
    start_block = run(r).events(strcmp(run(r).events.value, 'start_block'),:);
    gait_events = sortrows([gait_events; start_block]);
    % first remove invalid gait events:
    % ...remove gait_events (with 5-sec margin) that do not fall within the valid window
    idx = find(~([gait_events.onset]-5>run(r).info.validwindow(1) & [gait_events.onset+gait_events.duration]+5 < run(r).info.validwindow(2)));
    gait_events = gait_events(setdiff(1:end, idx),:);
    excl_gait = excl_gait + length(idx);
    % ...remove gait_events that fall within unexpected movements (5-sec
    % margin)
    if any(strcmp(run(r).events.type, 'Unexpected_movements'))
        unexp_mov = run(r).events(strcmp(run(r).events.type, 'Unexpected_movements'),:);
        for i=1:height(unexp_mov)
            idx = overlappingevt(gait_events, unexp_mov.onset(i)-5, unexp_mov.onset(i)+unexp_mov.duration(i)+5);
            gait_events = gait_events(setdiff(1:end, idx),:);
            excl_gait = excl_gait + length(idx);
        end
    end
    
    % loop over gait_events and collect info
    gait_eventsFT = table2FT(gait_events, sf);
    gait_trl = [];
    FOG_events = run(r).events(strcmp(run(r).events.type, 'FOG_Trigger'),:);
    for i=1:height(gait_events)
        % exclude gait_events with close FOG events within a 10 seconds margin (so if
        % FOG before, during or after normal gait event trial)
        idx = overlappingevt(FOG_events, gait_events.onset(i)-margin_FOG, gait_events.onset(i) + gait_events.duration(i) + margin_FOG);
        if ~isempty(idx)
            excl_gait = excl_gait + 1;
            continue
        else
            trlbegin = gait_eventsFT(i).sample - offset;
            trlend = gait_eventsFT(i).sample + offset;
            type = 2;
            [value, ~] = value_lookup(gait_eventsFT(i).value, gait_events(i,:));
            offset_trigger = nan;
            distance_trigger = nan;
            trial_duration = gait_events.duration(i);
            trl_new = [trlbegin trlend -offset type value offset_trigger distance_trigger trial_duration r];
            gait_trl = [gait_trl; trl_new];
        end
    end
    
    %%  calculate relevant movement data
    cfg = [];
    cfg.channel = find(contains(run(r).data_motion.data_raw.label, {'seg_Pelvis_velocity', 'seg_Pelvis_acceleration', 'jC1Head_jointAngle_'}));
    data_motion = ft_selectdata(cfg, run(r).data_motion.data_raw);
    
    % remove artifact of PD63
    if strcmp(ID, 'PD63') & r==4
        data_motion.trial{1}([1:3],11737:11763) = nan;
    end
    
    % convert pelvis velocity & acceleration from global to local reference frame
    cfg = [];
    cfg.channel = find(contains(run(r).data_motion.data_raw.label, 'seg_Pelvis_orientation'));
    data_rot = ft_selectdata(cfg, run(r).data_motion.data_raw);
    rotm = quat2rotm(data_rot.trial{1}');
    velo_new = nan(3, length(data_motion.time{1}));
    accel_new = nan(3, length(data_motion.time{1}));
    for i = 1:length(data_motion.time{1})
        velo_new(:,i) = rotm(:,:,i)'*data_motion.trial{1}(1:3,i);
        accel_new(:,i) = rotm(:,:,i)'*data_motion.trial{1}(4:6,i);
    end
    data_motion.trial{1}(1:3,:) = velo_new;
    data_motion.trial{1}(4:6,:) = accel_new;
    
    % low pass filter to remove step artefacts
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 0.7;
    data_motion_filt = ft_preprocessing(cfg, data_motion);
    
    % mirror the pitch (so negative = down) & roll (from left to right) of
    % the head orientation
    data_motion_filt.trial{1}(9,:) = -data_motion_filt.trial{1}(9,:);
    data_motion_filt.trial{1}(7,:) = -data_motion_filt.trial{1}(7,:);
    
    %% define trials
    cfg=[];
    cfg.trl=[FOG_trl; gait_trl];
    if isempty(cfg.trl)
        continue
    end
    if z_score
        tmp_nirs(r).data = ft_redefinetrial(cfg, run(r).data_nirs.data_zscore);
    else
        tmp_nirs(r).data = ft_redefinetrial(cfg, run(r).data_nirs.data_preproc);
    end
    % add movement data
    tmp_motion(r).data = ft_redefinetrial(cfg, data_motion_filt);
    
end % loop over runs

% boxplot of trigger_offsets
% if ~isempty(offset_trigger_all)
%   figure; boxplot(offset_trigger_all(:,2), offset_trigger_all(:,1)); hold on;
%   yline(0); title(ID)
% end


%% append data into one data_epoch
if isempty(tmp_nirs)
    data_epoch_nirs = [];
    data_epoch_motion = [];
else
    tmp_cmp_nirs=[tmp_nirs.data]; % get rid of empty structures
    data_epoch_nirs=tmp_cmp_nirs(1);
    tmp_cmp_motion=[tmp_motion.data]; % get rid of empty structures
    data_epoch_motion=tmp_cmp_motion(1);
    cfg=[];
    cfg.keepsampleinfo='yes';
    for i=2:length(tmp_cmp_nirs)
        data_epoch_nirs=ft_appenddata(cfg, data_epoch_nirs, tmp_cmp_nirs(i));
        data_epoch_motion=ft_appenddata(cfg, data_epoch_motion, tmp_cmp_motion(i));
    end
    if isfield(data_epoch_nirs, 'hdr')
        data_epoch_nirs = rmfield(data_epoch_nirs, 'hdr'); % otherwise conflict in ft_databrowser
    end
end

%% visualize
if ~isempty(data_epoch_nirs) & vis
    cfg                = [];
    cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
    cfg.viewmode       = 'butterfly';
    cfg.continuous = 'no';
    cfg.channel = {'Rx*'};
    cfg.fontsize=5;
    cfg.nirsscale =0.5;
    %   cfg.ylim = [-1 1];
    cfg.linecolor = 'rb';
    cfg.colorgroups = repmat([1 2],1, length(data_epoch_nirs.label)/2)';
    ft_databrowser(cfg, data_epoch_nirs);
    pause;
end

% show number of excluded events
warning('Excluded %d FOG events for this subject ', excl_FOG)
warning('Excluded %d gait events for this subject', excl_gait)


%% SUBFUNCTIONS
function [idx] = overlappingevt(annotations, begintime, endtime)
% find the indices of annotation events that fall within the event with the
% given [begintime endtime].

idx=find(([annotations.onset]<=begintime & [annotations.onset + annotations.duration]>begintime) |... % annotation includes the begintime
    ([annotations.onset]<endtime & [annotations.onset + annotations.duration]>=endtime) | ... % annotation includes the endtime
    ([annotations.onset]>=begintime & [annotations.onset + annotations.duration]<endtime)) ; % annotation falls within the event

function [value, trig_event] = value_lookup(FOG_value, trig_event)
if contains(FOG_value, '180')
    if contains(FOG_value, 'R')
        trig_event = trig_event(contains(trig_event.value, '180_R'),:);
        if contains(FOG_value, 'walking')
            value = 111;
        elseif contains(FOG_value, 'standing')
            value = 112;
        end
    elseif contains(FOG_value, 'L')
        trig_event = trig_event(contains(trig_event.value, '180_L'),:);
        if contains(FOG_value, 'walking')
            value = 121;
        elseif contains(FOG_value, 'standing')
            value = 122;
        end
    end
elseif contains(FOG_value, {'Doorway', 'doorway'})
    trig_event = trig_event(contains(trig_event.value, 'doorway'),:);
    if contains(FOG_value, 'walking')
        value = 21;
    elseif contains(FOG_value, 'standing')
        value = 22;
    end
elseif contains(FOG_value, {'Target', 'stop'})
    trig_event = trig_event(contains(trig_event.value, {'180', 'stop_walking', 'doorway'}),:);
    if contains(FOG_value, 'turn')
        value = 3.1;
    elseif contains(FOG_value, 'door')
        value = 3.2;
    else
        value = 3;
    end
elseif contains(FOG_value, {'SH', 'start_walking'})
    trig_event = trig_event(contains(trig_event.value, {'start_walking', '180'}),:);
    value = 4.1;
elseif contains(FOG_value, 'start_block')
    trig_event = 'NA';
    value = 4.2;
else
    value = 5;
end
if isempty(trig_event) & value~=5
    error('No relevant trigger event found for this FOG')
end

function [angle] = calculate_turnangle(data_motion, r, ID, trigger_event, FOG_event)
% calculates the turn angle at which FOG turn occurs

cfg = [];
cfg.channel = find(contains(data_motion.label, 'seg_Pelvis_orientation'));
ang_pelvis = ft_selectdata(cfg, data_motion);
% convert orientations to euler angles
[x, y, z]=q2e(ang_pelvis.trial{1}(1,:), ang_pelvis.trial{1}(2,:),ang_pelvis.trial{1}(3,:),ang_pelvis.trial{1}(4,:)); % convert orientation to euler angles
orient=rad2deg([x;y;z]);
orient_z = unwrap(orient(3,:), 340);
% define turn angle
angle_beginturn = orient_z(trigger_event.sample);
angle_endturn = orient_z(trigger_event.sample + trigger_event.duration);
if abs(angle_endturn-angle_beginturn) > 190 | abs(angle_endturn-angle_beginturn) < 160
    warning('no full 180 degree turn detected: %.0f', abs(angle_endturn-angle_beginturn))
    %         figure; hold on; plot(orient_z)
    %         plot(trigger_event.sample, angle_beginturn, 'ro')
    %         plot(trigger_event.sample + trigger_event.duration, angle_endturn, 'ro')
end
if FOG_event.sample < trigger_event.sample
    % FOG occured before start of turn; calculate the distance from the
    % turn instead
    angle = calculate_doordistance(data_motion, r, ID, trigger_event, FOG_event);
elseif FOG_event.sample > trigger_event.sample + trigger_event.duration
    angle = 1;
else
    angle = (orient_z(FOG_event.sample)-angle_beginturn)/(angle_endturn - angle_beginturn);
end

function distance_door = calculate_doordistance(data_motion, r, ID, trigger_event, FOG_event)
% Calculates the distance from the door at which the FOG doorway occurs
[rot_tra] = rotation_translation(data_motion, r, ID, 'visualize', false);
data_COM=reframe(data_motion, 'seg_COM_centerOfMass', 'rot+tra', rot_tra);
pos_door = data_COM.trial{1}(1, trigger_event.sample);
pos_FOG = data_COM.trial{1}(1, FOG_event.sample);
if FOG_event.sample<=trigger_event.sample
    distance_door = -abs(pos_door-pos_FOG);
else
    distance_door = abs(pos_door-pos_FOG);
end
