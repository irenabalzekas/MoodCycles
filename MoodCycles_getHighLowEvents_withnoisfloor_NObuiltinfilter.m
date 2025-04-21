% Function
% IB 12/15/2022

% Overall goal is to get the phase in the timeseries data of the
% observations in obs. Ie in what phase of the spike cycle do the low mood
% scores appear?

% INPUTS
% obs = values of observed events, ie mood scores
% obstime = timestamps of the observations
% threshold = percentage of points at the high and low end of all events to
%        consider, ie the top 25% of scores and the bottom 25% of mood
%        scores. To use the 25%, enter the number 4. 
% phase  = phase timeseries


% OUTPUTS
% highphases = phase of the top events
% lowphases = phase of the low events
% highfloor = randomized phase of the high end (noise floor). result of 100
%       repetitions. Reshuffle obs in time 100x, then re-sample high and
%       low. 
% lowfloor = randomized phase of the low end (noise floor). result of 100
%       repetitions
% ----- the following come from z-scored(obs) -------------
% lowsamps = values of the low samples
% lowtimes = timestamps of the low samples
% highsamps = values of the high samples
% hightimes = timestamps of the high samples



function [highphases, lowphases, allval_phases, random_highphases, random_lowphases, lowsamps, lowtimes, highsamps, hightimes] = MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter(phase, time, obs, obstime, threshold)

    % define low/high mood, above/below median
    quart = round(length(obs)/threshold); % number of samples in a quarter of the data
    sorted = sort(obs);
    lowthresh = sorted(quart); % 25th percentile and below
    highthresh = sorted(end-quart); % 75th percentile and above

    
    % find indices in mood time of low/high mood
    lowtimes = obstime(find(obs <= lowthresh));
    hightimes = obstime(find(obs >= highthresh));

    obs_z = zscore(obs);
    lowsamps = obs_z(find(obs <= lowthresh));
    highsamps = obs_z(find(obs >= highthresh));

    % get index in data of the low and high obs events
    lowind = [];
    for i = 1: length(lowtimes)
        [val,loc] = min(abs(time - lowtimes(i)));
        lowind(i) = loc;
    end 
    highind = [];
    for i = 1: length(hightimes)
        [val,loc] = min(abs(time - hightimes(i)));
        highind(i) = loc;
    end 
    % get index of all the obs events 
    obsind = [];
    for i = 1: length(obstime)
        [val,loc] = min(abs(time - obstime(i)));
        obsind(i) = loc;
    end 
    
    % to ensure overhang doesnt mean events at the ends are over sampled
    lowind = unique(lowind);
    highind = unique(highind);
    obsind = unique(obsind);
    

    % Get phase per high/low
    lowphases = phase(:,lowind);
    highphases = phase(:,highind);
    allval_phases = phase(:,obsind);

    % Reshuffling/finding noise floor
    random_lowphases = [];
    random_highphases = [];
    
    for f = 1
        
        randomlowphases_temp = [];
        randomhighphases_temp = [];
        
        for ii = 1:100 % repeating 100x

            random_obs = obs(randperm(length(obs))); % shuffle mood in time
            random_lowtimes = obstime(find(random_obs <= lowthresh));% get high and low indices from the shuffled
            random_hightimes = obstime(find(random_obs >= highthresh));

            % find indices in data vector of these events in time
            r_lowind = [];
            for i = 1: length(random_lowtimes)
                [val,loc] = min(abs(time - random_lowtimes(i)));
                r_lowind(i) = loc;
            end 
            r_highind = [];
            for i = 1: length(random_hightimes)
                [val,loc] = min(abs(time - random_hightimes(i)));
                r_highind(i) = loc;
            end 

            % Get phase per high/low
            randomlowphases_temp(ii,:) = phase(f,r_lowind);
            randomhighphases_temp(ii,:) = phase(f,r_highind);
            
        end
        
        random_lowphases(:,:,f) = randomlowphases_temp;
        random_highphases(:,:,f) = randomhighphases_temp;
    end 
   
end 
