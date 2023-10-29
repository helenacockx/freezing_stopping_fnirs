function run = nirs_preprocessing(run)
% NIRS_PREPROCESSING
% This function performs the preprocessing steps on the fNIRS data
%
% INPUT:
% - run: matlab structure containing the raw nirs & motion data of each run
%
% OUTPUT:
% - run: same as input, but with the preprocessed data in the field
% run.data_nirs.data_preprocessed
%
% DEPENDENCIES:
% - detect_MA.m
% - corr_MAspline.m
% - corr_MAwavelet.m
% - shortchannel_regression.m
% - auxchannel_regression.m

%%
method = 'ACR'; % perform auxilliary channel regression (so regression with short channels and motion data of the head)
MA_method = 'spline + wavelet';
vis = false;
plot_data = false;

time_MA = 0; time_total = 0;
for r=[1:length(run)]
    if isempty(run(r).data_nirs) | isempty(run(r).data_motion)
        continue
    else
        data_nirs = run(r).data_nirs.data_raw;
        data_motion = run(r).data_motion.data_raw;
    end
    fprintf('----- Processing data of run %d...  -----\n', r)
    
    %% calculate the step frequency and task frequency
    % step frequency
    idx_leftheel = find(strcmp(run(r).data_motion.data_raw.label, 'fc_LeftFoot_Heel_footContacts'));
    idx_rightheel = find(strcmp(run(r).data_motion.data_raw.label, 'fc_RightFoot_Heel_footContacts'));
    footsteps_L = ([diff(run(r).data_motion.data_raw.trial{1}(idx_leftheel,:)) 0])==1; % heel contacts
    footsteps_R = ([diff(run(r).data_motion.data_raw.trial{1}(idx_rightheel,:)) 0])==1; % heel contacts
    all_steps = (footsteps_L | footsteps_R);
    n = length(all_steps);
    Y = fft(all_steps);
    P2 = abs(Y/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = run(r).data_motion.data_raw.fsample*(0:n/2)/n;
    i_select = find(f<2.5);
    [M, I] = max(P1(10:i_select(end))); % find local max
    step_freq = f(I);
    
    % task frequency
    gait_events = run(r).events(strcmp(run(r).events.type, 'gait_event'),:);
    %     figure; plot(gait_events.onset, ones(size(gait_events.onset)), '+');
    task_freq = 1/mean(diff(gait_events.onset));
    task_freq_sd = std(1/diff(gait_events.onset));
    
    if vis
        figure('Name', sprintf('Freq. plot run %d', r)); plot(f(i_select), P1(i_select)); hold on; plot(f(I),M, 'o');
        xline(task_freq, 'g', 'LineWidth', 2);
        ax = gca;
        ft_plot_box([0.008 0.1 ax.YLim(1) ax.YLim(2)], 'tag', 'artifact', 'edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.1);
        ft_plot_box([task_freq-2*task_freq_sd task_freq+2*task_freq_sd ax.YLim(1) ax.YLim(2)], 'tag', 'artifact', 'edgecolor', 'none', 'facecolor', 'g', 'facealpha', 0.3);
        legend({'Power spectrum', 'step frequency', 'task frequency', 'band pass filter', 'task frequency 2 SD'})
        title(sprintf('Step frequency of %.02f Hz and task frequency of %.03f +- %.03f Hz in run %d \n', step_freq, task_freq, task_freq_sd, r));
    end
    if step_freq < 0.1
        warning('Step frequency lower than low pass filter')
    end
    if  task_freq > 0.5 | task_freq < 0.05
        warning('Task frequency close to band pass filter')
    end
    
    %% remove bad data segments and select motion data
    % remove bad data segments
    cfg = [];
    if strcmp(run(r).info.complete, 'no')
        cfg.latency = run(r).info.validwindow;
    end
    cfg.channel = data_nirs.label(find(~contains(data_nirs.label, 'ADC')));
    data_nirs_trim = ft_selectdata(cfg, data_nirs); % remove ADC channels
    
    % select motion data of the head
    cfg.channel = find(contains(data_motion.label, 'seg_Head_acceleration'));
    data_accel = ft_selectdata(cfg, data_motion);
    cfg.channel = find(contains(data_motion.label, 'seg_Head_angularAcceleration'));
    data_gyr = ft_selectdata(cfg, data_motion);
    cfg.channel = find(contains(data_motion.label, 'jC1Head_jointAngle_')); % jC1Head_jointAngle_
    data_orient = ft_selectdata(cfg, data_motion);
    % convert acceleration and angular acceleration from global to local reference frame
    cfg.channel = find(contains(data_motion.label, 'seg_Head_orientation'));
    data_rot = ft_selectdata(cfg, data_motion);
    rotm = quat2rotm(data_rot.trial{1}');
    accel_new = nan(3, length(data_accel.time{1}));
    gyr_new = nan(3, length(data_gyr.time{1}));
    for i = 1:length(data_accel.time{1})
        accel_new(:,i) = rotm(:,:,i)'*data_accel.trial{1}(:,i);
        gyr_new(:,i) = rotm(:,:,i)'*data_gyr.trial{1}(:,i);
    end
    data_accel_new = data_accel;
    data_accel_new.trial{1} = accel_new;
    data_gyr_new = data_gyr;
    data_gyr_new.trial{1} = gyr_new;
        
    %% motion artifact correction
    addpath('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\TDDR')
    if contains(MA_method, 'spline')
        [art] = detect_MA(data_nirs_trim, run(r).info.bad_channels);
        if ~isempty(art)
            % calculate % artifact
            time_MA = time_MA + sum(art(:,2)-art(:,1));
            time_total = time_total + length(data_nirs_trim.time{1})*length(data_nirs_trim.label);
            data_corr = corr_MAspline(data_nirs_trim, art);
        else
            data_corr = data_nirs_trim;
        end
        idx_good_chan = find(~contains(data_corr.label, run(r).info.bad_channels));
    else
        data_corr = data_nirs_trim;
    end
    if contains(MA_method, 'wavelet')
        data_corr = corr_MAwavelet(data_corr);
    end
    
    if vis
        figure; plot(data_corr.trial{1}(idx_good_chan,:)');
    end
    
    if r==length(run)-1
        warning('%.04f of the data was detected and corrected as motion artifacts', time_MA/time_total)
    end
    
    %% first band pass filter (of nirs and motion data)
    cfg =[];
    cfg.lpfilter = 'yes';
    cfg.lpfilterord = 2;
    cfg.lpfreq = 0.5;
    data_nirs_bp = ft_preprocessing(cfg, data_corr);
    data_accel_bp = ft_preprocessing(cfg, data_accel_new);
    data_gyr_bp = ft_preprocessing(cfg, data_gyr_new);
    data_orient_bp = ft_preprocessing(cfg, data_orient);
    
    %% OD to conc. changes   
    cfg                 = [];
    cfg.target          = {'O2Hb', 'HHb'};
    cfg.channel         = 'nirs';
    data_conc           = ft_nirs_transform_ODs(cfg, data_nirs_bp);
    
    %% Remove the bad channels 
    cfg = [];
    cfg.channel = data_conc.label(find(~contains(data_conc.label, run(r).info.bad_channels)));
    data_pruned = ft_selectdata(cfg, data_conc);
    cfg = [];
    cfg.latency = [0 data_conc.time{1}(end)];
    data_accel_pruned = ft_selectdata(cfg, data_accel_bp);
    data_gyr_pruned = ft_selectdata(cfg, data_gyr_bp);
    data_orient_pruned = ft_selectdata(cfg, data_orient_bp);
    
    %% systemic filtering
    switch method        
        % short channel regression
        case 'SCR'
            cfg = [];
            cfg.method = 'OLS';
            cfg.PCA = true; % first performs a PCA on the short channels and only uses the first 8 componenents for the regression
            data_filt = shortchannel_regression(cfg, data_pruned);
        % auxilliary channel regression with short channels and motion data
        case 'ACR'
            data_filt= auxchannel_regression(data_pruned, data_accel_pruned, data_gyr_pruned, data_orient_pruned);
        % no regression    
        case 'none'
            data_filt = data_pruned;
    end
    
    %% extra band pass filter after systemic filtering
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfiltord = 2;
    cfg.hpfreq = [0.01];
    cfg.lpfilter =  'yes';
    cfg.lpfiltord = 6;
    cfg.lpfreq = [0.1];
    data_filt_bpf = ft_preprocessing(cfg, data_filt);
    
    %% fill the bad channels with nans & save
    % make sure the channels are in the same order as the original data (just
    % to be sure that nothing get intermixed)
    data_new = data_conc;
    idx_good_chan = find(contains(data_conc.label, data_filt_bpf.label));
    data_new.trial{1}(idx_good_chan,:) = data_filt_bpf.trial{1}(:,:);
    idx_bad_chan = find(contains(data_conc.label, run(r).info.bad_channels));
    data_new.trial{1}(idx_bad_chan,:) = nan;
    % save under data_preproc
    run(r).data_nirs.data_preproc = data_new;
    
    %% visualize
    if plot_data
        for xi=1:length(data_conc.label)
            if ismember(xi, idx_good_chan)
                figure;
                ax(1)=subplot(3,1,1);hold on;
                title(data_conc.label{xi})
                plot(data_nirs_trim.time{1}, data_nirs_trim.trial{1}(xi,:)');
                plot(data_corr.time{1}, data_corr.trial{1}(xi,:)')
                plot(data_nirs_bp.time{1}, data_nirs_bp.trial{1}(xi,:)', 'Linewidth', 1)
                legend({'original', 'MA correction', 'low pass filter'});
                ax(2)=subplot(3,1,2);hold on;
                plot(data_conc.time{1}, data_conc.trial{1}(xi,:)');
                legend({'HbO'});
                ax(3)=subplot(3,1,3); hold on;
                idx_new = find(strcmp(data_conc.label{xi}, data_filt.label));
                plot(data_filt.time{1}, data_filt.trial{1}(idx_new, :)');
                plot(data_filt_bpf.time{1}, data_filt_bpf.trial{1}(idx_new,:)')
                legend({'ACR', 'bpf'})
                linkaxes(ax, 'x')
            end
        end
        figure; hold on;
        plot(data_accel_bp.time{1}, data_accel_bp.trial{1}'); plot(data_gyr_bp.time{1}, data_gyr_bp.trial{1}'), plot(data_orient_bp.time{1}, data_orient_bp.trial{1}')
        legend({'accel x', 'accel y', 'accel z', 'gyr x', 'gyr y', 'gyr z', 'orient x', 'orient y', 'orient z'});
        idx_SC = find(contains(data_conc.label, {'a ', 'b ', 'c ', 'd '}));
        [coeff,score,latent,tsquared,explained] = pca(data_conc.trial{1}(idx_SC,:)');
        figure; hold on;
        plot(data_conc.time{1}, score(:,1:8));
        legend();
    end
end


