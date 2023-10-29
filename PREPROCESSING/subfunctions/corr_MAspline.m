function [data_corr] = corr_MAspline(data, artifact)
% CORR_MASPLINE
% This function corrects the motion artifacts with a spline interpolation
% and is an adapted version of the Homer 3 function:
% hmrR_MotionCorrectSpline (p = 0.99)
%
% INPUT:
% - data: the fNIRS data in fieldtrip structure
% - artifact: artifact fieldtrip structure containing all the detected
% motion artifacts
%
% OUTPUT:
% - data_corr: the fNIRS data with spline interpolated data at the data
% segments of the motion artefacts
%
%%
vis = 0; % visualize output?

d = data.trial{1}';
dSpline = d;
t = data.time{1}';
fs = data.fsample;

p = 0.99;

% window widths limits for computing the mean in the segment shifts
dtShort = 0.3;  % seconds
dtLong  = 3;    % seconds

for ii = 1:length(data.label)
  lstMA = artifact(artifact(:,3)==ii,:); % sublist of motion artifact segments
  
  if ~isempty(lstMA)
    lstMs = lstMA(:,1); % starting indexes of mvt segments
    lstMf = lstMA(:,2);  % ending indexes of mvt segments
    lstMl = lstMf-lstMs;    % lengths of MA segments
    nbMA = length(lstMl);   % number of MA segments
    
    %% Do the spline interpolation on each MA segment
    for jj = 1:nbMA
      lst = lstMs(jj):(lstMf(jj)-1);
      % spline interp
      SplInterp = csaps(t(lst)', d(lst,ii)', p, t(lst)')';
      % corrected signal = original signal - spline interpolation
      dSpline(lst,ii) = d(lst,ii) - SplInterp;
    end
    
    %% Reconstruction of the whole time series (shift each segment)
    % First MA segment: shift to the previous noMA segment if it exists,
    % to the next noMA segment otherwise
    lst = (lstMs(1)):(lstMf(1)-1);
    SegCurrLength = lstMl(1);
    if SegCurrLength < dtShort*fs
      windCurr = SegCurrLength;
    elseif SegCurrLength < dtLong*fs
      windCurr = floor(dtShort*fs);
    else
      windCurr = floor(SegCurrLength/10);
    end
    
    if lstMs(1)>1
      SegPrevLength = length(1:(lstMs(1)-1));
      if SegPrevLength < dtShort*fs
        windPrev = SegPrevLength;
      elseif SegPrevLength < dtLong*fs
        windPrev = floor(dtShort*fs);
      else
        windPrev = floor(SegPrevLength/10);
      end
      meanPrev = mean(dSpline(lst(1)-windPrev:(lst(1)-1), ii));
      meanCurr = mean(dSpline(lst(1):(lst(1)+windCurr-1), ii));
      dSpline(lst,ii) = dSpline(lst,ii) - meanCurr + meanPrev;
      
    else
      if length(lstMs)>1
        SegNextLength = length(lstMf(1):(lstMs(2)));
      else
        SegNextLength = length(lstMf(1):size(d,1));
      end
      if SegNextLength < dtShort*fs
        windNext = SegNextLength;
      elseif SegNextLength < dtLong*fs
        windNext = floor(dtShort*fs);
      else
        windNext = floor(SegNextLength/10);
      end
      meanCurr = mean(dSpline((lst(end)-windCurr):(lst(end)-1),  ii));
      meanNext = mean(dSpline((lst(end)+1):(lst(end)+windNext), ii));
      dSpline(lst,ii) = dSpline(lst,ii) - meanCurr + meanNext;
    end
    
    
    % Intermediate segments
    for kk=1:(nbMA-1)
      % no motion
      lst = lstMf(kk):(lstMs(kk+1)-1);
      SegPrevLength = lstMl(kk);
      SegCurrLength = length(lst);
      if SegPrevLength < dtShort*fs
        windPrev = SegPrevLength;
      elseif SegPrevLength < dtLong*fs
        windPrev = floor(dtShort*fs);
      else
        windPrev = floor(SegPrevLength/10);
      end
      if SegCurrLength < dtShort*fs
        windCurr = SegCurrLength;
      elseif SegCurrLength < dtLong*fs
        windCurr = floor(dtShort*fs);
      else
        windCurr = floor(SegCurrLength/10);
      end
      meanPrev = mean(dSpline((lst(1)-windPrev):(lst(1)-1), ii));
      meanCurr = mean(d(lst(1):(lst(1)+windCurr-1), ii));
      
      dSpline(lst,ii) = d(lst,ii) - meanCurr + meanPrev;
      
      % motion
      lst = (lstMs(kk+1)):(lstMf(kk+1)-1);
      SegPrevLength = SegCurrLength;
      SegCurrLength = lstMl(kk+1);
      if SegPrevLength < dtShort*fs
        windPrev = SegPrevLength;
      elseif SegPrevLength < dtLong*fs
        windPrev = floor(dtShort*fs);
      else
        windPrev = floor(SegPrevLength/10);
      end
      if SegCurrLength < dtShort*fs
        windCurr = SegCurrLength;
      elseif SegCurrLength < dtLong*fs
        windCurr = floor(dtShort*fs);
      else
        windCurr = floor(SegCurrLength/10);
      end
      meanPrev = mean(dSpline((lst(1)-windPrev):(lst(1)-1), ii));
      meanCurr = mean(dSpline(lst(1):(lst(1)+windCurr-1), ii));
      
      dSpline(lst,ii) = dSpline(lst,ii) - meanCurr + meanPrev;
    end
    
    % Last not MA segment
    if lstMf(end)<length(t)
      lst = (lstMf(end)-1):length(t);
      SegPrevLength = lstMl(end);
      SegCurrLength = length(lst);
      if SegPrevLength < dtShort*fs
        windPrev = SegPrevLength;
      elseif SegPrevLength < dtLong*fs
        windPrev = floor(dtShort*fs);
      else
        windPrev = floor(SegPrevLength/10);
      end
      if SegCurrLength < dtShort*fs
        windCurr = SegCurrLength;
      elseif SegCurrLength < dtLong*fs
        windCurr = floor(dtShort*fs);
      else
        windCurr = floor(SegCurrLength/10);
      end
      meanPrev = mean(dSpline((lst(1)-windPrev):(lst(1)-1), ii));
      meanCurr = mean(d(lst(1):(lst(1)+windCurr-1), ii));
      
      dSpline(lst,ii) = d(lst,ii) - meanCurr + meanPrev;
    end
    
    if vis
      figure; hold on; plot(d(:, ii)); plot(dSpline(:, ii));
      ax = gca;
      lim_y = [ax.YLim(1)-0.1 ax.YLim(2)+0.1];
      for jj = 1:nbMA
        ft_plot_box([lstMs(jj) lstMf(jj) lim_y], 'tag', 'artifact', 'edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.1);
      end
      ylim(lim_y);
      title(data.label{ii});
    end
  end
end

data_corr = data;
data_corr.trial{1} = dSpline';
    
    
    