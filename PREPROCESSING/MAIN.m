%% MAIN
% Preprocesses the nirs data and saves this for each participant.
%
% INPUT:
% - *run.mat files (contains raw nirs & motion data)
% - *events.mat files (contains the events)
% - participants.xlsx: table with participant demographical and clinical
% info
%
% OUTPUT:
% - *run.mat files (now also containing the processed data:
% run.data_nirs.data_preproc and/or run.data_nirs.data_zscore)
% - *trails*.mat files: nirs data (*trialsnirs_z.mat) and motion data
% (*trialsmotion.mat)segmented in trials based on the events
%
% DEPENDENCIES:
% - load_run_info.m; load_events.m; remove_badchannels.m;
% nirs_preprocessing.m; define_trials.m
% - age2DPF.m
% - Homer3 toolbox


%% script to analyze the nirs data
root_dir='C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data';

% add path to scripts
addpath(genpath(fullfile(root_dir, 'scripts','final')))
addpath(genpath('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\Homer3\FuncRegistry\UserFunctions'))
addpath('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\fieldtrip\external\artinis')
dbstop if error

% run following processing steps
proc_info = 0;
proc_events = 0;
proc_dpf = 0;
proc_badchan = 0;
proc_preproc = 0;
proc_zscore = 0;
proc_trial = 1;
use_zscore = 0;

%% run over participant groups
group = {'PD', 'HC'};
for g = 1:length(group)
    switch group{g}
        case 'PD'
            ID={'PD06', 'PD10', 'PD11', 'PD12', 'PD15', 'PD16', 'PD17', 'PD22', 'PD25',...
                'PD31', 'PD35', 'PD45', 'PD46', 'PD48', 'PD50', 'PD57', 'PD61', 'PD62', 'PD63',...
                'PD76', 'PD77', 'PD88',  'PD90', 'PD96'};
        case 'HC'
            ID={'HC01', 'HC02', 'HC04', 'HC06', 'HC13', 'HC19', 'HC21', 'HC28', 'HC29', 'HC33', 'HC34', 'HC35'...
                'HC41', 'HC42', 'HC57', 'HC60', 'HC66', 'HC67', 'HC68', 'HC76', 'HC90', 'HC91'};
    end
    ID=sort(ID);
    
    % read in participants info
    participants = readtable(fullfile(root_dir, 'processed/participants.xlsx'), 'Sheet', 'Summary');
    
    %% run over participants
    for i=1:length(ID)
        % load data & events
        fprintf('====== sub-%s ======= \n', ID{i})
        load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_run.mat', ID{i})));
        load(fullfile(root_dir, 'processed', ['sub-' ID{i}], sprintf('sub-%s_events.mat', ID{i})));
        
        %% load run info
        if proc_info
            [run] = load_run_info(run, ID{i});
        end
        
        %% load events from ELAN
        if proc_events
            [run]=load_events(run, events, ID{i});
        end
        
        %% update DPF info based on participant age and wavelength
        if proc_dpf
            % load participant information
            idx = find(strcmp(participants.ID, ID{i}));
            for r=1:length(run)
                if ~isempty(run(r).data_nirs)
                    age = participants.age(idx);
                    for c=1:length(run(r).data_nirs.data_raw.opto.label)
                        wl_idx = nonzeros(run(r).data_nirs.data_raw.opto.tra(c,:));
                        if wl_idx(1) == -wl_idx(2)
                            wl = run(r).data_nirs.data_raw.opto.wavelength(abs(wl_idx(1)));
                        else
                            error
                        end
                        DPF = age2DPF(age, wl);
                        run(r).data_nirs.data_raw.opto.DPF(c) = DPF;
                    end
                end
            end
        end
        
        %% detect bad channels
        if proc_badchan
            run = remove_badchannels(run, ID{i});
        end
        %
        %% nirs preprocessing
        if proc_preproc
            run = nirs_preprocessing(run);
        end
        
        %% z-score
        if proc_zscore
            for r=1:length(run)
                if ~isempty(run(r).data_nirs) & ~isempty(run(r).data_motion)
                    run(r).data_nirs.data_zscore = run(r).data_nirs.data_preproc;
                    run(r).data_nirs.data_zscore.trial{1} = zscore(run(r).data_nirs.data_preproc.trial{1}, 0, 2);
                end
            end
        end
        %% trial definition
        if proc_trial
            [trials_nirs, trials_motion] = define_trials(run, ID{i}, 0, use_zscore);
        end
        
        %% save run & events
        mkdir(fullfile(root_dir, 'processed', 'final',['sub-' ID{i}]))
        if use_zscore
            save(fullfile(root_dir, 'processed', 'final',['sub-' ID{i}], sprintf('sub-%s_trialsnirs_z.mat', ID{i})), 'trials_nirs');
        else
            save(fullfile(root_dir, 'processed', 'final',['sub-' ID{i}], sprintf('sub-%s_trialsnirs.mat', ID{i})), 'trials_nirs');
        end
        save(fullfile(root_dir, 'processed', 'final', ['sub-' ID{i}], sprintf('sub-%s_trialsmotion.mat', ID{i})), 'trials_motion');
    end
end

