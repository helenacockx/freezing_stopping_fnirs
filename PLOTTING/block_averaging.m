%% BLOCK_AVERAGING
% Performs a first quick analysis based on the grand average of the nirs data
% over participants (without Bayesian statistics). 
% 
% INPUT: 
% - sub-*_trialsnirs.mat: the epoched fNIRS data of each participant
%
% OUTPUT:
% - plots (and videos) of the block average of each participant and plots of the grand
% average over participants
%
% DEPENDENCIES:
% - timelock.m 
% - grandaverage.m 

%% generate general across-subject block averaging of nirs data
clear data_timelock
group = {'PD', 'HC'};
type = {'normal'}; % options: FOG or normal
condition = {'turn walking', 'doorway walking', 'start', 'stop'};

for g=1:length(group)
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
    
    for t=1:length(type)
        for c=1:length(condition)
            for i=1:length(ID)
                load(fullfile(root_dir, 'processed', 'final', ['sub-' ID{i}], sprintf('sub-%s_trialsnirs.mat', ID{i})));
                data_timelock{i} = timelock(trials_nirs, ID{i}, 0, type{t}, condition{c}, true);
            end
            grandaverage(data_timelock, type{t}, condition{c}, group{g});
        end
    end
end

