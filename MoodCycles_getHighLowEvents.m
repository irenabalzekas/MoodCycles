% INPUTS
% obs = values of observed events, ie mood scores
% obstime = timestamps of the observations
% threshold = percentage of points at the high and low end of all events to
%        consider, ie the top 25% of scores and the bottom 25% of mood
%        scores. To use the 25%, enter the number 4. 

% OUTPUTS
% lowtimes = timestamps of the low samples
% hightimes = timestamps of the high samples


function [lowind, highind] = MoodCycles_getHighLowEvents(threshold,obs,obstime)

    % define low/high mood, above/below median
    quart = round(length(obs)/threshold); % number of samples in a quarter of the data
    sorted = sort(obs);
    lowthresh = sorted(quart); % 25th percentile and below
    highthresh = sorted(end-quart); % 75th percentile and above

    lowind = find(obs <= lowthresh);
    highind = find(obs >= highthresh);
    
%     % find indices in mood time of low/high mood
%     lowtimes = obstime(find(obs <= lowthresh));
%     hightimes = obstime(find(obs >= highthresh));
% 
%     lowsamps = obs(find(obs <= lowthresh));
%     highsamps = obs(find(obs >= highthresh));
    

end