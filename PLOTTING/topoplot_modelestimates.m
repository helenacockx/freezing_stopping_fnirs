%% TOPOPLOT_MODELESTIMATES
% This script creates topoplots based on the model estimates of each
% channel.
%
% INPUT:
% - summary/*/*_channels.csv: summary output table from R containing model
% estimates of each channel for the given condition
%
% OUTPUT:
% - topoplots of the model estimates for each of the conditions

%% Initialisation
fig_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\figures\final\fNIRS\model_estimates';
figure_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\paper\figures\topoplots';
conditions = {'turn', 'door', 'start', 'stop', 'standing', 'walking'};

% load layout
load('C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\freezing_fnirs\data\scripts\layout.mat');
layout.width = 2*layout.width;
layout.height = 2*layout.height;

%% create a topoplot for each condition
for c=1:length(conditions)
    % read summary table
    if contains(conditions{c}, 'FOG')
        file = ['\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\scripts\final\summary\FOG\' conditions{c} '_channel.csv'];
        factor = {'Intercept'};
    else
        file = ['\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\scripts\final\summary\nrl\' conditions{c} '_channel.csv'];
        factor = {'Intercept','group1', 'HC', 'PD'};
    end
    t = readtable(file);
    
    for f=1:length(factor)
        % create FT structure
        idx = find(strcmp(t.factor, factor{f}));
        data.avg = t.Estimate(idx);
        data.label = t.channel(idx);
        data.time = [0];
        prob = t.Prob(idx);
        if strcmp(factor{f}, 'group1')
            data.avg = -data.avg;
        end
        
        % create topoplot
        cfg =[];
        cfg.parameter = 'avg';
        cfg.xlim = [0 3];
        cfg.zlim = [-0.3 0.3];
        cfg.layout = layout;
        cfg.style = 'straight';
        cfg.interplimits = 'mask_individual';
        cfg.interpolation = 'v4';
        cfg.shading = 'interp';
        cfg.gridscale = 300;
        cfg.colorbar = 'no';
        cfg.comment = 'no';
        cfg.highlight = repmat({'on'}, length(data.label), 1);
        cfg.highlight(prob<0.975) = {'off'};
        cfg.highlightsize = num2cell((prob-0.95)*350);
        cfg.highlightchannel = cellfun(@(x) {x}, data.label, 'UniformOutput', false);
        cfg.highlightsymbol = repmat({'p'}, length(data.label), 1);
        cfg.highlightcolor = num2cell((data.avg<0)*[1 1 1], 2);
        ft_topoplotER(cfg, data)
        % title(factor{f});
        % mkdir(fullfile(fig_dir, conditions{c}))
        % saveas(gcf, fullfile(fig_dir, conditions{c}, sprintf('%s.jpg', factor{f})))
        
        % save plot for figure in paper
        cfg.highlightsize = num2cell((prob-0.95)*250);
        ft_topoplotER(cfg, data)
        axis tight
        set(gca, 'fontsize', 10)
        set(gcf, 'Units', 'centimeters')
        set(gcf,'Position',[1 1 9 9]);
        mkdir(fullfile(figure_dir, conditions{c}))
        filename = fullfile(figure_dir, conditions{c}, factor{f});
%         print(filename, '-depsc', '-r600')
        
    end
end