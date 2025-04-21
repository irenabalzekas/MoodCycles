% Function for least squares linear phase filter design
% IB 2/13/22, implementation from Nick Gregg

% INPUT
% cycleranges = list of cycle windows ie  1-5 days, 20-30 days. Each window
        % will be converted into a filter family with 20 total frequencies within
        % that range 
% signal = timeseries to filter
% fs = sampling rate. samples/second of signal

% OUTPUT
% filtered = firls noncausal filtered signal for each band 
% phase = hilbert of each filtered


function [filtered, phase] = firls_nick(signal, fs, cycleranges)
    
    filtered = zeros(length(cycleranges), length(signal));
    phase = zeros(length(cycleranges), length(signal));

    for c = 1: length(cycleranges)
        filtfam = linspace(cycleranges{c}(1), cycleranges{c}(2),20); % convert cycle ranges to filter families around central frequency - scale density based on cycle length/central frequency
        
        % INITIALIZE AND DEFINE FILTER PARAMETERS
        f = nan(length(filtfam), 6);
        % convert fs to nyquist framework (samples per day/2)
        nyq = round(fs *(60*60*24)/2);
        n = nyq*3;  %order-n FIR (how many points in filter) (generally btw 2 and 5, larger order better freq discrimintation, not as tight temporally)
        for i = 1:length(filtfam)
            f(i,:) = cat(2, [.025, .5, .75, 1.25, 1.5]/nyq/filtfam(i), 1);% f = [.025, .5, .75, 1.25, 1.5, nyq]/nyq; %MikeCohenPg182. freq characteristics of filter relative to nyquist freq.. (baseline; start transition up; end transition; end plateau; transition down; nyquist freq)
        end 
        A = [0, 0, 1, 1, 0, 0];  %stop and pass components of filter
 
        % MAKE FILTER KERNEL
        for i = 1:length(filtfam)
            filterFIR(i).multi=firls(round(n*filtfam(i)), f(i,:), A);
        end 
        
        % FILTER DATA (non-causal (filtfilt)) 
        data_firls = zeros(length(filtfam), length(signal));
        for j = 1:length(filtfam) 
            data_firls(j,:) = filtfilt(filterFIR(j).multi, 1, signal)*sqrt(filtfam(j)); 
        end
        
        % ADD UP ALL THE FILTERED TIMESERIES WITHIN THE CYCLE RANGE FOR THE FINAL OUTPUT
        %filtered(c,:) = sum(data_firls,1);
        filtered(c,:) = mean(data_firls,1);

        phase(c,:) = angle(hilbert(sum(data_firls,1)));
    end 

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Nick's original code 
% %%  filter params
% % firls FIR filter.  field 1: #samples. field 2: [
% filtfam = cat(2, [.3:.1:2],[2.5:.5:10], [11:1:40]); % this was center frequencies for the filter family in days
% 
%  
% %% initialize and define some parameters
% f = nan(length(filtfam), 6);
% nyq = 12; %for hourly samples nyq = 12; for q6hr samples nyg = 2; q3hr nyq=4;
% n = nyq*3;  %order-n FIR (how many points in filter) (generally btw 2 and 5, larger order better freq discrimintation, not as tight temporally)
% % f = [.025, .5, .75, 1.25, 1.5, nyq]/nyq; %MikeCohenPg182. freq characteristics of filter relative to nyquist freq.. (baseline; start transition up; end transition; end plateau; transition down; nyquist freq)
% for i = 1:length(filtfam)
%     f(i,:) = cat(2, [.025, .5, .75, 1.25, 1.5]/nyq/filtfam(i), 1);
% end 
% A = [0, 0, 1, 1, 0, 0];  %stop and pass components of filter
%  
% % % % % filterFIR1=firls(n, f, A);
% % % % % data_firls1 = filtfilt(filterFIR1, 1, data);
%  
% %% allocate
% data_avgfiltBig = zeros(length(filtfam),length(data));
% data_firlsmulti = zeros(length(filtfam),length(data));
% data_firlsmulti_causal = zeros(length(filtfam),length(data));
% x_multi = zeros(length(filtfam),length(data));
% x_multicausal = zeros(length(filtfam),length(data));
% x_runBig = zeros(length(filtfam),length(data));
%  
% %% make filter kernel
% for i = 1:length(filtfam)
%     filterFIR(i).multi=firls(round(n*filtfam(i)), f(i,:), A);
% end 
%  
% %% filter data; non-causal (filtfilt) and causal (filter) & hilbert analytic signal
% % for j = 1:length(filtfam) %OG
% for j = 1:3%length(filtfam)
% 
%     data_firlsmulti(j,:) = filtfilt(filterFIR(j).multi, 1, data)*sqrt(filtfam(j));
%     data_firlsmulti_causal(j,:) = filter(filterFIR(j).multi, 1, data)*sqrt(filtfam(j));
% %     data_firlsmulti(j,:) = filtfilt(filterFIR(j).multi, 1, data(:,2))*sqrt(filtfam(j));
% %     data_firlsmulti_causal(j,:) = filter(filterFIR(j).multi, 1, data(:,2))*sqrt(filtfam(j));
%     x_multi(j,:) = hilbert(data_firlsmulti(j,:));
%     x_multicausal(j,:) = hilbert(data_firlsmulti_causal(j,:));
%     x_runBig(j,:) = hilbert(data_avgfiltBig(j,:)); 
% end
%  
% %% power and phase data from hilbert analytic signal
% phi_runBig = angle(x_runBig);
% % phi_circ = angle(x_circ);
% phi_multi = angle(x_multi);
% phi_multicausal = angle(x_multicausal);
%  
% power_runBig = abs(x_runBig);
% % power_circ = abs(x_circ);
% power_multi = abs(x_multi);
% power_multicausal = abs(x_multicausal);filter
%  
% %% Sz phase locking - find the phase at which each seizure occurs for each frequency in the filter family
% SzT_IEAtime = zeros(1,numel(SzT));
%  
% for i = 1:numel(SzT)
%     [SzT_IEAtime(i), idx(i)] = min(abs(data(:,3)-SzT(i)));
% end
% 
