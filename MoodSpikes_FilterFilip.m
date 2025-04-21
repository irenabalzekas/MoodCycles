% Function to bandpass filter in the really low frequency range (cycles around 100 days)
% Filter designed by Filip Mivalt, December 2022

% INPUT
% x = signal vector
% t = time vector (timestamps, in seconds)
% rng1 = [high_end, low_end] vector of the cycle range **in units of days**
% order = filter order. defaults to 10001. % higher orders e.g. 50001 will provide better filtering. However will take longer. 10001 should work as well

% OUTPUT
% x_filtered = filtered signal

function [x_filtered] = MoodSpikes_FilterFilip(x, t,rng1, order)

    
%     t = spikeinfo.time;
%     x = spikeinfo.signal;
    fs = 1 / (mean(diff(t)) / 3600);
    pxx = abs(fft(x)) / length(x);
    pxx = pxx(1:length(x)/2).^2;
    f = linspace(0, fs/2, length(pxx));
%     rng1 = spikeinfo.cycleranges{2}; % pick range that you want
    min_cutoff = min(rng1); % extending the band +- 10 might help;
    max_cutoff = max(rng1); 
    % % min_cutoff = 60;
    % % max_cutoff = 110; % custom range
    min_cutoff = 1 / (min_cutoff * 24); % convert to frequency in hr^-1
    max_cutoff = 1 / (max_cutoff * 24);

%     % plot frequency response 
%     figure(1)
%     loglog(f, pxx)
%     xline(min_cutoff, 'r')
%     xline(max_cutoff, 'r')
%     xlabel('Frequency hr^-^1')
%     ylabel('pow -')

    if isempty(order) == 0
        order = order;
    elseif isempty(order) == 1
        order = 10001;
    end 

    b_hp = fir1(order, max_cutoff*(fs/2)); % get filters
    b_lp = fir1(order, min_cutoff*(fs/2));
    a = 1;
 
    [pxx_hp, w] = freqz(b_hp, a, 2*length(f)); % estimate the transition characteristic
    pxx_hp = abs(pxx_hp);
    b_hp = b_hp / max(pxx_hp); % change the gain of the filter to 1 exactly in time domain
    pxx_hp = pxx_hp / max(pxx_hp); % in frequency domain
    [pxx_lp, w] = freqz(b_lp, a, 2*length(f)); % same for this one
    pxx_lp = abs(pxx_lp);
    b_lp = b_lp / max(pxx_lp);
    pxx_lp = pxx_lp / max(pxx_lp);
    pxx_final = pxx_lp.*(1-pxx_hp); % estimate the final bandpass
% 
%     % plot the spectra
%     figure(2)
%     loglog(w/2/pi*fs, 1-pxx_hp) % plot the spectra
%     hold on;
%     loglog(w/2/pi*fs, pxx_lp)
%     loglog(w/2/pi*fs, pxx_final)
%     legend('HP', 'LP', 'BP')
%     hold off;

    zero_pend = zeros(1, 2*order); % zeros that will be appended to the signal to proceed with the filtration
    x_pend = [zero_pend, x, zero_pend];
    x_pend = filtfilt(b_lp, a, x_pend); % filter LP
    x_pend = x_pend - filtfilt(b_hp, a, x_pend); % filter HP
    x_filtered = x_pend(length(zero_pend)+1:end-length(zero_pend)); % remove the appended zeros - might be 1 sample off
    x_filtered = x_filtered / max(pxx_final); % correct for the gain of the bandpass created by combination of both filters

%     % plot original and filtered 
%     figure(3);
%     plot(x)
%     hold on;
%     plot(x_filtered)
%     hold off; 

end 


%% Original demo from Filip Mivalt, December 2022
% 
% subject = 'M4';
% filename = strcat('P:\Personal\Irena\MoodCyclesSpikeCycles_Paper\realdata\',subject,'\spikes\spikeinfo.mat');
% load(filename)
% 
% t = spikeinfo.time;
% x = spikeinfo.signal;
% fs = 1 / (mean(diff(t)) / 3600);
% pxx = abs(fft(x)) / length(x);
% pxx = pxx(1:length(x)/2).^2;
% f = linspace(0, fs/2, length(pxx));
% rng1 = spikeinfo.cycleranges{2}; % pick range that you want
% min_cutoff = min(rng1); % extending the band +- 10 might help;
% max_cutoff = max(rng1); 
% % % min_cutoff = 60;
% % % max_cutoff = 110; % custom range
% min_cutoff = 1 / (min_cutoff * 24); % convert to frequency in hr^-1
% max_cutoff = 1 / (max_cutoff * 24);
% 
% % plot frequency response 
% figure(1)
% loglog(f, pxx)
% xline(min_cutoff, 'r')
% xline(max_cutoff, 'r')
% xlabel('Frequency hr^-^1')
% ylabel('pow -')
% 
% order = 10001; % higher orders e.g. 50001 will provide better filtering. However will take longer. 10001 should work as well
% b_hp = fir1(order, max_cutoff*(fs/2)); % get filters
% b_lp = fir1(order, min_cutoff*(fs/2));
% a = 1;
%  
% [pxx_hp, w] = freqz(b_hp, a, 2*length(f)); % estimate the transition characteristic
% pxx_hp = abs(pxx_hp);
% b_hp = b_hp / max(pxx_hp); % change the gain of the filter to 1 exactly in time domain
% pxx_hp = pxx_hp / max(pxx_hp); % in frequency domain
% [pxx_lp, w] = freqz(b_lp, a, 2*length(f)); % same for this one
% pxx_lp = abs(pxx_lp);
% b_lp = b_lp / max(pxx_lp);
% pxx_lp = pxx_lp / max(pxx_lp);
% pxx_final = pxx_lp.*(1-pxx_hp); % estimate the final bandpass
% 
% % plot the spectra
% figure(2)
% loglog(w/2/pi*fs, 1-pxx_hp) % plot the spectra
% hold on;
% loglog(w/2/pi*fs, pxx_lp)
% loglog(w/2/pi*fs, pxx_final)
% legend('HP', 'LP', 'BP')
% hold off;
% 
% zero_pend = zeros(1, 2*order); % zeros that will be appended to the signal to proceed with the filtration
% x_pend = [zero_pend, x, zero_pend];
% x_pend = filtfilt(b_lp, a, x_pend); % filter LP
% x_pend = x_pend - filtfilt(b_hp, a, x_pend); % filter HP
% x_filtered = x_pend(length(zero_pend)+1:end-length(zero_pend)); % remove the appended zeros - might be 1 sample off
% x_filtered = x_filtered / max(pxx_final); % correct for the gain of the bandpass created by combination of both filters
% 
% % plot original and filtered 
% figure(3);
% plot(x)
% hold on;
% plot(x_filtered)
% hold off; 