  function eventFT = table2FT(eventT, fsample)
%   converts event in table format with timestamps to events in Fieldtrip
%   format (structure format) with samples based on the given sample
%   frequency
  eventT.Properties.VariableNames={'timestamp', 'duration', 'type', 'value'};
  eventT.duration=round(eventT.duration*fsample); % convert to samples
  eventT.sample=round(eventT.timestamp*fsample+1);
  eventFT=table2struct(eventT);
  