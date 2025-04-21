% generate a sinusoidal signal with chosen sinusoids, noise, polynomial trend, 

% OUTPUT
% years = number of years signal duration
% samprate = samples per hours
% maxsigamp = scaling factor for component sinusoids
% poi = periods of interest (period, in days, of component sinusoids)
% polycoef = polynomial coefficients 
% noiseparams = added random noise scaling factor for each component
% shiftup = 1 (if set ==1, then values are shifted to all be greater than
% zero (ie can't have negative spike rate)

% OUTPUT
% simsig = signal
% t = time, timestamp (seconds) per sample
% n = number of samples in signal
% dt = time between samples 


function [simsig, t, n, dt] = simsignal(years, samprate, maxsigamp, poi, polycoef, noiseparams, shiftup)

    n = samprate * 24 * 365 * years; % samples per hour for one year 
    totalsec = years * 60 * 60 * 365 * 24;
    dt = round(totalsec/n); % year in seconds/n - % confirm this is indeed 20 min***

    % SIMULATE SIGNAL WITH NOISE
    Fs = 1/dt; % sampling period and rate are inverse 
    f = Fs*(0:n-1)/n;
    f_adj = 1./(f*24*60*60); 
    foi = [];
    for i = 1: length(poi)
        [val, idx] = min(abs(f_adj-poi(i)));
        foi(i) = f(idx);
    end 

    t = linspace(1, totalsec, n);
    basesig = zeros(1, n);
    components = [];
    noisecomponents = [];
    for i= 1:length(foi)
        component = maxsigamp(i).*sin(2*pi*foi(i)*t);
        addednoise = noiseparams(i).*randn(1,n);
        basesig = basesig + component + addednoise;
        components(i,:) = component;
        noisecomponents(i,:) = addednoise;
    end
    t4pol = [1:n];
    polynomialcomponent = polycoef(1)*t4pol + polycoef(2)*t4pol.^2 + polycoef(3)*t4pol.^3 +polycoef(4)*t4pol.^4 +polycoef(5)*t4pol.^5;
    simsig = basesig + polynomialcomponent; % adding noise
    
    if shiftup == 1
        simsig = simsig + abs(min(simsig)); % shift up so all values above zero and keep some around zero
    else
        simsig = simsig;
    end 
end

