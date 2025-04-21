% IB 12/15/22

% INPUT
% events = timestamps of an event of interest (ie seizures)
% phase = timeseries of phase associated with overall timeseries of
% interest

% OUTPUT
% eventphases = phase in data timeseries of the events in the event
%   timeseries - filtered according to cycle range


function [eventphases] = MoodCycles_getEventPhase_Nobuiltinfilter(events,time, phase)

    szind = [];
    for i = 1: length(events)
        [val,loc] = min(abs(time - events(i)));
        szind(i) = loc;
    end 
    szind = unique(szind); 

    eventphases = phase(:,szind);

end 