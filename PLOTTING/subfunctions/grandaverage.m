function grandaverage(data_timelock, type, condition, group)
%% GRANDAVERAGE
% Computes the block average over multiple participants for a given
% condition and visualizes this in multiple plots and videos.
% 
% INPUT: 
% - data_timelock: cell-array with timelocked data of all participants
% (output from timelock.m)
% - type: type of the condition of the given timelocked data
% ('FOG' or 'normal')
% - condition: condition of the given timelocked data
% (e.g., 'turn R walking', 'stop turn', 'doorway')
% - group: group of participants for who to calculate the grand average
% ('PD' or 'HC')
%
% OUTPUT:
% - plots and videos.
%
% DEPENDENCIES:
% - plotCI.m

%% Initialisation
video = 0; % to create a video of the topoplots over time (1/0)
CI_plots = 0; % to create plots of HbO and HbR of each channel with confidence intervals (1/0)

figure_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\figures\final\fNIRS';
% figpaper_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\paper\figures';

% load layout
load('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\scripts\layout.mat');
layout.width = 2*layout.width;
layout.height = 2.5*layout.height;

%% Calculations
% remove empty subjects
emptysubj=find(cellfun(@isempty,data_timelock));
subjects=setdiff(1:length(data_timelock), emptysubj);

% remove the DPF fields, otherwise ft_timelockgrandaverage starts to
% complain that the opto information is not unique
for s=1:length(subjects)
  if isfield(data_timelock{subjects(s)}.opto, 'DPF')
    data_timelock{subjects(s)}.opto = rmfield(data_timelock{subjects(s)}.opto, 'DPF');
  end
end

% calculate the grand average
cfg=[];
data_GA=ft_timelockgrandaverage(cfg,data_timelock{subjects});

%% Plot
% first separate O2Hb and HHb channels
cfg=[];
cfg.channel='* [O2Hb]';
data_O2Hb=ft_selectdata(cfg, data_GA);
% and rename labels such that they have the same name as HHb channels
for i=1:length(data_O2Hb.label)
  tmp = strsplit(data_O2Hb.label{i});
  data_O2Hb.label{i}=tmp{1};
  data_O2Hb.opto.label{i} = tmp{1};
end
cfg=[];
cfg.channel='* [HHb]';
data_HHb=ft_selectdata(cfg, data_GA);
% and rename labels such that they have the same name as HHb channels
for i=1:length(data_HHb.label)
  tmp = strsplit(data_HHb.label{i});
  data_HHb.label{i}=tmp{1};
  data_HHb.opto.label{i} = tmp{1};
end
% then plot both on the lay-out
cfg                   = [];
cfg.channel           = find(~contains(data_O2Hb.label, {'a', 'b', 'c', 'd'}));
cfg.showlabels        = 'no';
cfg.layout            = layout;
cfg.showoutline       = 'yes';
cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor        = 'rb'; % O2Hb is showed in red, HHb in blue
% cfg.linecolor = [0.83529  0.36863  0];
cfg.linewidth       = 2;
cfg.ylim              = [-0.05 0.1];
if strcmp(type, 'FOG')
  cfg.ylim = [-0.1 0.1]; % alternatives: [-0.1 0.1]; [-0.2 0.3]; [-0.5 0.5]
end
cfg.xlim              = [-6 10];
cfg.skipcomnt  = 'yes';
cfg.skipscale = 'yes';
ft_multiplotER(cfg, data_O2Hb, data_HHb); 
title([type ' ' condition ' ' group])
set(gcf, 'PaperOrientation', 'portrait');
filename = fullfile(figure_dir, 'timecourse', sprintf('%s_%s_%s.jpg', type, condition, group))
saveas(gcf, filename);

% figure for the paper
cfg.linewidth       = 4;
ft_multiplotER(cfg, data_O2Hb, data_HHb); 
set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [1 1 20 20]);
% filename = fullfile(figpaper_dir, 'timecourses', sprintf('%s_%s_%s', type, condition, group));
% saveas(gcf, [filename '.png'])

% % plot short channel data
% cfg = [];
% cfg.channel = {'*a*', '*b*', '*c*', '*d*'};
% SC_data_O2Hb{k} = ft_selectdata(cfg, data_O2Hb{k});
% SC_data_HHb{k} = ft_selectdata(cfg, data_HHb{k});
% cfg = [];
% cfg.layout            = layout;
% cfg.showlabels = 'no';
% cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
% cfg.linecolor        = 'rb'; % O2Hb is showed in red, HHb in blue
% cfg.linewidth       = 2;
% cfg.ylim              = [-1 1];
% cfg.xlim              = [-6 10];
% ft_multiplotER(cfg, SC_data_O2Hb{k}, SC_data_HHb{k}); 
% % ft_multiplotER(cfg, SC_data_O2Hb{k}); 

% topoplot
cfg =[];
cfg.channel = {'all', '-*a', '-*b', '-*c', '-*d'};
cfg.parameter = 'avg';
cfg.xlim = [0 3];
cfg.zlim = [-0.05 0.05]; % alternative: [-0.05 0.05]; [-0.3 0.3]
cfg.layout = layout;
cfg.style = 'straight';
cfg.interplimits = 'mask_individual';
cfg.interpolation = 'v4';
cfg.shading = 'interp';
cfg.gridscale = 100;
% cfg.colorbar = 'yes';
cfg.comment = 'no';
ft_topoplotER(cfg, data_O2Hb)
title([type ' ' condition ' ' group])
filename = fullfile(figure_dir, 'topoplots', sprintf('topoplotHbO_%s_%s_%s.jpg', type, condition, group));
saveas(gcf, filename);
ft_topoplotER(cfg, data_HHb)
title([type ' ' condition ' ' group])
filename = fullfile(figure_dir, 'topoplots', sprintf('topoplotHbR_%s_%s_%s.jpg', type, condition, group));
saveas(gcf, filename);

%% create video 
if video
  cd(fullfile(figure_dir, 'videos'));
  begtime = -10;
  increment = 1; % stepwise in seconds
  
  % make a figure with the desired size, see https://en.wikipedia.org/wiki/Display_resolution
  figh = figure('Color',[1 1 1]);
  set(figh, 'WindowState', 'maximized');
%   title([type ' ' condition ' ' group])
  
  % prepare the video file
  vidName = sprintf('topoplotHbO_%s_%s_%s', type, condition, group);
  vidObj = VideoWriter(vidName, 'MPEG-4');
  vidObj.FrameRate = 1/increment;
  vidObj.Quality = 100;
  open(vidObj);
  
  while begtime<=9
    cfg.xlim = [begtime begtime+1];
    cfg.zlim = [-0.05 0.05];
    cfg.figure = figh;
    ft_topoplotER(cfg, data_O2Hb);
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
    % go to the next segment of data
    begtime = begtime + increment;
  end
end

%% plot O2Hb and HHb of each channel with confidence intervals
if CI_plots
  for i=1:length(data_O2Hb.label)
    figure;
    b=plotCI(data_HHb, data_HHb.label{i}, length(data_GA.cfg.previous), 'b');
    a=plotCI(data_O2Hb, data_O2Hb.label{i}, length(data_GA.cfg.previous), 'r');
    xline(0, '--', 'LineWidth', 2);
    xlabel('time (s)')
    legend([a b],{'HbO2', 'HHb'})
    title([type ' ' condition ' ' group ' ' data_O2Hb.label{i}])
  end
end
