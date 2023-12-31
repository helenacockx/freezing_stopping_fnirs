function data_rcs=shortchannel_regression(cfg, datain);
% SHORTCHANNEL_REGRESSION
% shortchannelregression performs reference channel subtraction for fNIRS
% data. It is an adaptation from ft_nirs_referencechannelsubstraction 
% (http://www.fieldtriptoolbox.org/reference/ft_nirs_referencechannelsubtraction/)
%
% Use as
%   outdata = shortchannel_regression(cfg, indata)
% where indata is nirs data and cfg is a configuration structure that
% should contain:
%   cfg.method        = string, 'regstat2', 'QR' or 'OLS' (default = 'QR')
%   cfg.verbose       = boolean, whether text output is desired (default = false)
%   cfg.nearest       = only use the closest short channel for the
%   regression (in this case, takes the short channel with the same
%   receiver as the long channel)
%

%%
% get the options
cfg.method        = ft_getopt(cfg, 'method', 'QR');
cfg.verbose       = ft_getopt(cfg, 'verbose', false);
cfg.nearest       = ft_getopt(cfg, 'nearest', false);
cfg.PCA           = ft_getopt(cfg, 'PCA', false);

% find short and long channel indexes
SC=find(contains(datain.label, {'a ', 'b ', 'c ', 'd '})& all(~isnan(datain.trial{1}(:,:)),2)); %index of all short channels that does not contain nans
LC=find(~contains(datain.label, {'a ', 'b ', 'c ', 'd '})); %index of all long channels

data_rcs			 = datain;
data_rcs.label = datain.label(LC);
for tr=1:numel(datain.trial)
  shallow		= datain.trial{tr}(SC,:);
  shallow		= bsxfun(@minus,shallow,mean(shallow,2)); % mean detrend
  shallowlabel = datain.label(SC);
  
  deep		= datain.trial{tr}(LC,:);
  deep		= bsxfun(@minus,deep,mean(deep,2)); % mean detrend
  deeplabel = datain.label(LC);

  time		= datain.time{tr};
  
  % Reference channel subtraction
  ndeep		= size(deep,1);
  signal		= NaN(size(deep));
  if ~cfg.nearest
    x			= shallow'; % use all short channels
    idx_SC=1:length(shallowlabel);
    if cfg.PCA
      [coeff,score,latent,tsquared,explained] = pca(x);
      display(sum(explained(1:8)))
      display(explained(1:8)')
      x = score(:,1:8);
    end
  end
  for dpIdx	= 1:ndeep
    y				= deep(dpIdx,:)';
    if all(isnan(y))
      continue
    end
    
    if cfg.nearest
        % find corresponding short channel
        deeplabel_split=strsplit(deeplabel{dpIdx}, {'Tx', ' '});
        idx_SC=find(startsWith(shallowlabel, deeplabel_split{1})) ;       
        x = shallow(idx_SC,:)';
    end
    
    switch (cfg.method)
      case 'regstat2'
        b				= regstats2(y,x,'linear',{'beta','r'});
        beta    = b.beta;
        res			= b.r;
        
      case 'QR'
        [Q,R] = qr(x,0);
        beta = R\(Q'*y);
        yhat = x*beta;
        res = y - yhat;
        
      case 'OLS'
        x2 = [repmat(1, size(x, 1), 1) x];
        beta = x2\y;
        yhat = x2*beta;
        res  = y - yhat;
        cfg.verbose;
        beta(1) = [];
        
      otherwise % it should never come here as we use ft_checkopt
        error('unrecognized method');
    end
    
    signal(dpIdx,:) = res';
    
    % sanity check of results
    if cfg.verbose & ~cfg.PCA
      fprintf('Found the following meaningful shallow channels for deep channel %s:', deeplabel{dpIdx});
        shIdx = find(beta>0.5);
        for s=1:numel(shIdx)
          fprintf('\n\t%s', shallowlabel{idx_SC(shIdx(s))})
        end
      fprintf('\n\n')
    end
  end
  
  % overwrite
  data_rcs.time{tr}	= time;
  data_rcs.trial{tr}	= signal;
end