% Matlab code accompanying Manuscript: “A feasibility study of multiday mood 
% cycles in three people with temporal lobe epilepsy undergoing ambulatory 
% multimodal neuromonitoring and deep brain stimulation”
% Author: Irena Balzekas, 2025
% Please contact balzekasirena@gmail.com with questions

% DEPENDENCIES + FUNCTIONS
% BPWP package https://github.com/irenabalzekas/BPWP
% CVX matlab package http://cvxr.com/cvx/
% Matlab wavelet toolbox https://www.mathworks.com/products/wavelet.html
% Matlab circstat toolbox https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
% Matlab Watsons U2 https://www.mathworks.com/matlabcentral/fileexchange/43543-pierremegevand-watsons_u2

% Main script and additional functions
% - MoodCyclesSeizureCyclesForPublication
% - simsignal() 
% - firls_nick()
% - MoodCycles_getEventPhase_Nobuiltinfilter()
% - MoodCycles_getHighLowEvents()
% - MoodSpikes_FilterFilip
% - MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter

%% Define file paths
dest = '';

%% Define some variables (for data simulation)
subject = 'S1';
years = 1; % total data duration for simulated data
timeshift = 1.7e9; % Arbitrary number of seconds to shift unix time 

%% Simulate and plot spike rate timeseries - one year total, seizures 
% Modified from previously published approach to simulate spike rate timeseries
% Balzekas et al. 
% https://github.com/irenabalzekas/BPWP/blob/main/GenSimSpikeData_script.m

% DEFINE BASIC PARAMETERS + PATHS
timezone = 'America/Chicago';
samprate = 3; % number of samples per hour
n = samprate * 24 * 365 * years; % samples per hour for one year 
totalsec = years * 60 * 60 * 365 * 24;
dt = round(totalsec/n); % year in seconds/n - % confirm this is indeed 20 min***
maxsigamp = [60, 20, 10, 10, 40, 25, 20]; % helps scale amplitude of compoent sinusoid (periods of interest)
poi = [1, 7, 15, 21, 30, 50, 100]; %  DEFINE PERIODS IN SIM SIGNAL periods (in days) of interest to include in simulated signal
polycoef = [-.005 1*10^-50 10^-100  0 0]; % for first, second, third, fourth, fifth, 
noiseparams = [16, 10, 10, 10, 6, 4, 4]; % MAKE THIS MATCH WITH NUMBER OF POI***

% Generate signal, do CWT decomp
[simsig, t, n, dt] = simsignal(years, samprate, maxsigamp, poi, polycoef, noiseparams,1);
t = t + timeshift; % Shift time to arbitrary unix relevant window
[meanpowerCWTL2norm, perioddays_cwt, WT, fwt] = wavelet_decomp_L2norm(simsig', 1/dt); 
meanpowerCWTreg = nanmean(abs(WT).^2 ,2);

% Simulate seizure samples - can be more playful here, for simulation, just
% randomly sampled 100 time samples from the spike timeseries

szind = randsample(length(t), 100);
sz = t(szind);

% Create spike struct
spikedata.spikes = simsig;
spikedata.time = t;
spikedata.dates = datetime(t, 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
spikedata.perioddays_cwt = perioddays_cwt;
spikedata.cwtpower = meanpowerCWTreg;
spikedata.timezone = timezone;
spikedata.sz = sz;

% Plot, review timeseries 
figure
subplot(1,4,[1,2,3])
plot(spikedata.dates, spikedata.spikes, 'Color', 'black')
xlabel('Time')
ylabel('Spike rate')
title('Simulated spike rate timeseries')
subplot(1,4,4)
plot(perioddays_cwt, meanpowerCWTreg)
xlabel('Period (days)')
ylabel('Power')
title('CWT')
xlim([-10 max(perioddays_cwt)])

savename = strcat(dest,subject, '_spikedata.mat');
save(savename, 'spikedata')

%% EMA data: Simultate timeseries, define variables, import/format

% ema.mat  - table of EMA scores with rows = entries, columns = subject,
% questionnaire name, datetimecompleted, utctimecompleted, column for each
% response to each item, then total score

timezone = 'America/Chicago';

% generate simulated EMA data 

% DEFINE BASIC PARAMETERS - can keep some, such as sampling rate similar to
% spike rate section, for convenience
maxsigamp = [2, 10, 6, 8, 3, 4]; % 
poi = [0.7, 1, 14, 15, 32, 110]; %  DEFINE PERIODS IN SIM SIGNAL periods (in days) of interest to include in simulated signal
polycoef = [-.0002 1*10^-10 0 0 0]; % for first, second, third, fourth, fifth, 
noiseparams = [5, 2, 2, 2, 3, 2]; % MAKE THIS MATCH WITH NUMBER OF POI***

% SIMULATE SIGNAL WITH NOISE (can use same function as to generate spike data, just pick different parameters. We will also downsample it)
[simsig, t, n, dt] = simsignal(years, samprate, maxsigamp, poi, polycoef, noiseparams,0);

% Convert t to arbitrary unix period
t = t + timeshift;

% Randomly sample 120 values 
emasamp = sort(randsample([1:length(t)],120));

% CREATE STRUCT
emadata.signal = simsig(emasamp); % change name to EMA/ema struct.... 
emadata.time = t(emasamp);
emadata.time_days = datetime(t(emasamp), 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
emadata.timezone = timezone;

% PLOT GENERATED TIMESERIES 
figure
plot(emadata.time_days, emadata.signal, '.', 'MarkerSize', 7)
xlabel('Time (days)')
ylabel('Spike rate')
title('Simulated EMA timeseries')

savename = strcat(dest,subject, '_emadata.mat');
save(savename, 'emadata')

%% Run BPWP with delta parameter sweep 10% cross validation 

% DEFINE FIXED PARAMS PER SUBJECT
numrepeats = 10; % number of times to recalculate MSE
percent2drop = 10; % for cross validation 90/10
maxdegree = 3;
maxperiod = 120;
minperiod = 10; 
sampletype = 'densesubrange';
subjects = {'S1'};
Ns = [1500,2000,2500];
deltas = 0:100:10000; % 0:50000:1500000; - adjust range based on parameter sweep output, looking for depict global minimum

for p = 1:length(subjects)  
    
    subject = subjects{p};
    
    for n = 1: length(Ns)
    
        N = Ns(n);

        signal = emadata.signal;
        time = emadata.time; 
        dt = round((max(time) - min(time))) / (N-1);
        t_full = 1:N;

        % If times between samples are less than resulting dt, that's a
        % challenge for this data. you end up with the same phi for multiple
        % samples. THEREFORE. some of those samples will have to be omitted 
        tocut = find(diff(time) <= dt);
        signal(tocut) = [];
        time(tocut) = [];

        tic
        % Building the dct basis can be slow if its large, good to do outside loop
        [f_k, desiredperioddays] = frequency_sampling(N,maxperiod, minperiod, dt, sampletype);
        [dct_basis, scalefactor] = DCT2_basis(N, f_k); % making the basis takes a long time. want to pre-define it before running BPDN for loops and loops 
        [poly_basis] = polynomial_basis(N, maxdegree);
        toc
        fprintf('basis is built')

        MSEpersamp = [];
        Xspersamp = [];
        zspersamp = [];
        measureddatapersamp = [];
        tpersamp = [];
        phipersamp = [];
        training_md_persamp = [];
        training_t_persamp = [];
        training_phi_persamp = []; 
        testing_md_persamp = [];
        testing_t_persamp = [];
        testing_phi_persamp = [];

        for s = 1

            % Define numsamp needed 
            numsamp = floor(length(signal)*(100-percent2drop)/100); % this should probably be what is there minus 10% 

            [measureddata_init, t_init, phi_init, training_md, training_t, training_phi, testing_md, testing_t, testing_phi] = BPDN_samplerealdata_traintest(signal,time,numsamp, dt, percent2drop, numrepeats);
            % debug fig
            % figure;plot(phi_init, measureddata_init);hold on;plot(training_phi(1,:), training_md(1,:),'*');hold on;plot(training_phi(5,:), training_md(5,:),'o')

            % Initialize variables 
            MSE = zeros(length(deltas), numrepeats);
            xs = zeros(length(deltas),numrepeats,N);
            zs = zeros(length(deltas),numrepeats, maxdegree+1);

            parfor dd = 1: length(deltas)

                delta = deltas(dd);
                fprintf('delta sweeps with n delta =')
                fprintf(num2str(length(deltas)))
                fprintf('currently on iteration')
                fprintf(num2str(dd))

                for nn = 1: numrepeats

                    % Training data
                    measureddata = training_md(nn,:);
                    t = training_t(nn,:);
                    phi = training_phi(nn,:);

                    % Testing data
                    measured_missed = testing_md(nn,:);
                    t_missed = testing_t(nn,:);
                    phi_missed = testing_phi(nn,:);

                     % Run BPDN
                    [x, z, ~, ~, ~, ~] = BPDN(delta, measureddata', phi, dct_basis, poly_basis);
                    [reconsig] = BPDN_reconsig(f_k, x,scalefactor, maxdegree, z, t_full);
                    xs(dd,nn,:) = x;
                    zs(dd,nn,:) = z;

                    % Get missed samples from reconstruction
                    measured_missed_recon = reconsig(phi_missed); 

                    % Calculate MSE
                    MSE(dd,nn) = immse(measured_missed, measured_missed_recon); 

                end

            end 

            MSEpersamp{s} = MSE;
            Xspersamp{s} = xs;
            zspersamp{s} = zs;
            measureddatapersamp{s} = measureddata_init;
            tpersamp{s} = t_init;
            phipersamp{s} = phi_init;
            training_md_persamp{s} = training_md;
            training_t_persamp{s} = training_t;
            training_phi_persamp{s} = training_phi; 
            testing_md_persamp{s} = testing_md;
            testing_t_persamp{s} = testing_t;
            testing_phi_persamp{s} = testing_phi;

        end
        
        deltaselection.subject = subject;
        deltaselection.signal = signal;
        deltaselection.time = time;
        deltaselection.deltas = deltas;
        deltaselection.numsamp = numsamp;
        deltaselection.note = 'Repeated cross validation for each combo of delta and sample density (number of samples). Each cell corresponds to cross validation MSE vs delta for a particular density';
        deltaselection.MSE = MSEpersamp;
        deltaselection.N = N;
        deltaselection.dt = dt;
        deltaselection.measuredata_init = measureddatapersamp;
        deltaselection.t_init = tpersamp;
        deltaselection.phi_init = phipersamp;
        deltaselection.training_md_persamp = training_md_persamp;
        deltaselection.training_md_persamp= training_md_persamp;
        deltaselection.training_phi_persamp = training_phi_persamp; 
        deltaselection.testing_md_persamp = testing_md_persamp;
        deltaselection.testing_t_persamp = testing_t_persamp;
        deltaselection.testing_phi_persamp = testing_phi_persamp;
        deltaselection.dct_basis = dct_basis;
        deltaselection.maxdegree = maxdegree;
        deltaselection.minperiod = minperiod;
        deltaselection.maxperiod = maxperiod;
        deltaselection.sampletype = sampletype;
        deltaselection.poly_basis = poly_basis;
        deltaselection.f_k = f_k;
        deltaselection.desiredperioddays = desiredperioddays;
        deltaselection.x = Xspersamp;
        deltaselection.z = zspersamp;
        deltaselection.scalefactor = scalefactor;
        deltaselection.t_full = t_full;
        deltaselection.percent2drop = percent2drop;

        deltaselection.moodrecontime = (deltaselection.t_full.*deltaselection.dt) + deltaselection.time(1);
        recondates = datetime(deltaselection.moodrecontime, 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
        recondates.TimeZone = emadata.timezone;
        deltaselection.recondates = recondates;
        moodsampledates = datetime(deltaselection.time, 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
        moodsampledates.TimeZone = emadata.timezone;
        deltaselection.moodsampledates = moodsampledates;
        
        savename = strcat(dest,subject,'MSE_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_traintestdata.mat');
        save(savename,'deltaselection');
    end
end 

%% Plot delta param sweep outputs - MSE vs delta value for different N

subjects = {'S1'};
colors = {'red','black','blue','cyan','green','magenta'};
fsize = 15;
percent2drop = 10;
numrepeats = 10;
fig = figure('Position',[10 10 1500 500])
Ns = [1500,2000,2500]; 

for i = 1:length(subjects)
    
    subplot(1,length(subjects),i)
    subject = subjects{i};

    for nn = 1: length(Ns)
        
        filename = strcat(dest,subject,'MSE_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(Ns(nn)),'_traintestdata.mat');    
        load(filename);
    
        x = deltaselection.deltas;
        dat = deltaselection.MSE{1}'; % Number of ‘Experiments’ In Data Set
        yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
        ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
        CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

        %plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
        CIs = yCI95+yMean;
        curve1 = CIs(1,:);
        curve2 = CIs(2,:);
        xs = [x, fliplr(x)];
        inbetween = [curve1,fliplr(curve2)];
        fill(xs,inbetween, colors{nn},'FaceAlpha',0.5,'DisplayName','','EdgeColor','none','HandleVisibility','off')
        hold on
        labelname = strcat('N = ',{' '},num2str(Ns(nn)));
        plot(x, yMean, 'Color',colors{nn}, 'LineWidth',2,'DisplayName',labelname{:}) % Plot Mean Of All Experiments
        hold on
        [row,col] = find(yMean == min(yMean));
        plot(x(col),yMean(col),'o','LineWidth',2,'DisplayName',strcat('d=',num2str(x(col))))
        hold on
        grid
    end
    
       
    xlabel('Delta')
    ylabel('MSE')
    title(subjects{i})

    legend
    set(gca,'FontSize',fsize)
    
end 

%% Update deltaselection structure with delta value that minimized MSE
% Min value output in previous figure, manually input here 

subject = 'S1';
percent2drop = 10;
N = 2000; % input optimal N for participant based on parameter sweep
numrepeats = 10;
savename = strcat(dest,subject,'MSE_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_traintestdata.mat');
load(savename)

deltaselection.delta = 800; % input here based on best performing N and delta from above figure. minimum mse per curve

save(savename,'deltaselection');

%% Get spike cycles from equivalent time period as mood - calculate CWT and save
% note - need to hard code for right hemisphere. 

subjects = {'S1'};

for i = 1:length(subjects)
    
    subject = subjects{i};
    
    % LOAD SUBJECTS MOOD DATA (for temporal reference)
    load(strcat(dest,subject,'_emadata.mat'));
    moodtime = emadata.time;
    
    % LOAD SUBJECT'S SPIKE DATA
    load(strcat(dest,subject,'_spikedata.mat'));

    % this spike source is from interpolated spikes (all uniform dt)
    signal = spikedata.spikes;
    time = spikedata.time;
    
    % FIND indices in spike series that correspond with mood series, cut spikes down to match
    moodstart = min(moodtime);
    moodend = max(moodtime); 
    [vals, spikestart] = min(abs(time - moodstart)); 
    [vale, spikeend] = min(abs(time - moodend));

    % Cut time down to specified duration
    signal = signal(spikestart:spikeend);
    time = time(spikestart:spikeend);
    
    dt_ogdata = min(diff(time));
    
    % Calculate classic wavelet 
    [meanpowerCWTL2norm,perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(signal, 1/dt_ogdata);
    meanpowerCWTreg = nanmean(abs(WT).^2,2);
    reconsigCWT = icwt(WT);
    [perioddays_cwt,meanpowerCWTreg, peakvalues, peaklocs, peakperiods] = BPDN_cwt2findcycles(signal, dt_ogdata); % redundant, but just for peak finding
    
    % Save content into struct 
    spikeinfo.subject = subject;
    spikeinfo.hemisphere = [];
    spikeinfo.time = time;
    spikeinfo.signal = signal;
    spikeinfo.sz = spikedata.sz;
    spikeinfo.timezone = spikedata.timezone;
    spikeinfo.meanpowerCWTL2norm = meanpowerCWTL2norm;
    spikeinfo.perioddays_cwt = perioddays_cwt;
    spikeinfo.WT = WT;
    spikeinfo.fwt = fwt;
    spikeinfo.meanpowerCWTreg = meanpowerCWTreg;
    spikeinfo.reconsigCWT = reconsigCWT;
    spikeinfo.peakvalues = peakvalues;
    spikeinfo.peaklocs = peaklocs;
    spikeinfo.peakperiods = peakperiods;
    
    spikeinfo.spikedates = datetime(spikeinfo.time, 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
    spikeinfo.szdates = datetime(spikeinfo.sz, 'ConvertFrom','posixtime','TimeZone','UTC','Format','yy,mm hh:mm');
    
    savename = strcat(dest, subject, '_spikeinfo.mat');
    save(savename,'spikeinfo');
end 

%% Run significant peaks and add to deltaselection struct based on optimal delta 

subjects = {'S1'};
Ns = [2000];
percent2drop = 10;
numrepeats = 10;

for i = 1:length(subjects)
    
    subject = subjects{i};
    
    N = Ns(i);
    
    % Load reference data from parameter sweep (N, dt, data, etc.)
    destfolder = dest; 
    load(strcat(dest,subject,'MSE_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_traintestdata.mat'));

    % Define basics
    delta = deltaselection.delta;
    maxdegree = deltaselection.maxdegree;
    sampletype = deltaselection.sampletype;
    signal = deltaselection.signal;
    time = deltaselection.time; 
    dt_ogdata = time(2)-time(1);
    dt = deltaselection.dt;
    N = deltaselection.N;
    t_full = deltaselection.t_full;
    f_k = deltaselection.f_k;
    desiredperioddays = deltaselection.desiredperioddays;
    dct_basis = deltaselection.dct_basis;
    scalefactor = deltaselection.scalefactor;
    poly_basis = deltaselection.poly_basis;
    % Select data based on desired num per day
    ind = 1;

    measureddata_wds = deltaselection.measuredata_init{ind}(:);
    time_wds =  deltaselection.t_init{ind}(:);
    phi = round((time_wds-min(time))/dt) + 1; % indices for when each time data were sampled

    % First, BPDN
    [x_wds, xreshuffled, xpercentiles, z_wds, zreshuffled, zpercentiles, reshuff_optval] = BPDN_wReshuffling_short(delta, measureddata_wds, phi, dct_basis, poly_basis);
    sighits = find(xpercentiles > 99);
    [reconsig_wds] = BPDN_reconsig(f_k, x_wds, scalefactor, maxdegree, z_wds, t_full);

    % Rescale power by frequency (% Change from l1 to l2 norm for cwt.  l2 norm reduces amplitude of high freq data, while l1 does not. (scale is inversely related to frequency))
    x2 = [];
    for i=1:length(x_wds) 
        x2(i) = x_wds(i).*1/sqrt(f_k(i));  
    end
    x_power_rescale = abs(x2);

    % Nick CWT - with L2 adjustment
    [meanpowerCWTL2norm,perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(signal, 1/dt_ogdata);
    reconsigCWT = icwt(WT);
    [perioddays_cwt,meanpowerCWTreg,pks,locs] = BPDN_cwt2findcycles(signal, dt_ogdata); % redundant, but just for peak finding

    % Modify time scale for plots
    time_full_days = (time - min(time))./(60*60*24);
    t_days = (time_wds - min(time_wds))./ (60*60*24);
    t_recon_days = (t_full*dt) ./ (60*60*24);

    % Saving data for fig
    figs.x = x_wds;
    figs.xreshuffled = xreshuffled;
    figs.xpercentiles = xpercentiles;
    figs.z = z_wds;
    figs.zpercentiles = zpercentiles;
    figs.reshuff_optval = reshuff_optval;
    figs.measureddata = measureddata_wds;
    figs.time = time_wds;
    figs.sighits = sighits;
    figs.reconsig = reconsig_wds;
    figs.x_power_rescale = x_power_rescale;
    figs.meanpowerCWTreg = meanpowerCWTreg;
    figs.perioddays_cwt = perioddays_cwt;
    figs.WT = WT;
    figs.fwt = fwt;
    figs.reconsigCWT = reconsigCWT;
    figs.pks = pks;
    figs.locs = locs;
    figs.time_full_days = time_full_days;
    figs.t_days = t_days;
    figs.t_recon_days = t_recon_days;
    deltaselection.fig = figs;

    % Append and save
    savename = strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat');
    save(savename,'deltaselection');
    fprintf('subject done')
end

%% Notes on evaluating significant hits 

% 1. DCT can yield negative coefficients. 
% 2. The method may yield some very small DCT coefficients, on the order of 4e-06,
% which still remain outliers when compared with the distribution from 
% the non-parametric reshuffling occurs. 
% Plotted as power (the square), the values are so small they show up as zero. 
%       The question of whether or not to omit/ignore such model outputs is likely
%       model and application specific. 

%% get mean and median, variance IMS scores, min/max gaps between scores, stim durations, scores per stim paradigm per subject
% Simulated file of stimulation parameter start and end dates not included
% here "stimfile"

clear
close all

subs = {'S1'};
Ns = [2000];

means = [];
numentries = [];
medians = [];
variances = [];
maxgapindays = [];
mediangapindays = [];
sdgapsindays = [];
totEMAscoreduringhighstim = [];% number of entries 
totEMAscoreduringlowstim = [];
medEMAscoreduringhighstim = [];
medEMAscoreduringlowstim = [];
meanEMAscoreduringhighstim = [];
meanEMAscoreduringlowstim = [];
wilcoxonP_EMAscoreLowvHighstim = [];
totaldayslowstim = [];
totaldayshighstim = [];
for s = 1: length(subs)
   
    subject = subs{s};
    percent2drop = 10;
    numrepeats = 10;
    N = Ns(s);
    percentilethreshold = 99;%

    % Load reference data from parameter sweep (N, dt, data, etc.) + sighits
    load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
    % Load stimulation parameters 
    tbase = min(deltaselection.time);% centralize time
    stimfile = ''; % 
    stim = readtable(stimfile);
    stimstarts = (stim.start_uutc- tbase) ./ (24*60*60);
    stimends = (stim.stop_uutc - tbase) ./ (24*60*60);
    
    % general measures, median, variance EMA Scores etc. 
    numentries(s) = length(deltaselection.signal);
    means(s) = mean(deltaselection.signal);
    medians(s) = median(deltaselection.signal);
    variances(s) = var(deltaselection.signal);
    
    maxgapindays(s) = max(diff(deltaselection.time))/(60*60*24); % longest gap between entries
    mediangapindays(s) = median(diff(deltaselection.time)/(60*60*24));
    sdgapsindays(s) = std(diff(deltaselection.time)/(60*60*24)); % SD of gaps
    
    % ensure only considering stim params during EMA sampling period - cut down stim data file
    stimmod = stim;
    stimmod.StimStartpos = posixtime(stimmod.StimStart);
    stimmod.StimEndpos = posixtime(stimmod.StimEnd);
 
    indstart = [];
    indend = [];
    inside = [];
    % find EMAs in the stim file.
    for a = 1: length(deltaselection.time)
        if (deltaselection.time(a) < posixtime(stimmod.StimEnd(end))) & (deltaselection.time(a) > posixtime(stimmod.StimStart(1))) 
            inside(a) = 1;
        else
            inside(a) = 0;
        end 
    end   
    st = deltaselection.time(min(find(inside ==1)));
    et = deltaselection.time(max(find(inside==1)));
    
    replacedstartind = [];
    for aa = 1: height(stimmod)
        if (st < posixtime(stimmod.StimEnd(aa))) & st > posixtime(stimmod.StimStart(aa))% find the first EMA score within the start of stim
            stimmod.StimStartpos(aa) = st;
            replacedstartind = aa;
        else
            continue
        end 
    end 
    replacedendind = [];
    for aa = 1: height(stimmod)
        if (et < posixtime(stimmod.StimEnd(aa))) & (et > posixtime(stimmod.StimStart(aa))) % find the first EMA score within the start of stim
            stimmod.StimEndpos(aa) = et;
            replacedendind = aa;
        else
            continue
        end 
    end 
    %subselect just rows from stimmod beyond replacedstartind and replacedindend
    stimmod = stimmod(replacedstartind:replacedendind,:); % might be wrong
  
    % find total duration of time stim <50 Hz (low freq) and >= 50 Hz (high freq)
    lowtime = [];
    hightime = [];
    for i = 1: height(stimmod)
        if stimmod.Frequency(i) < 50
            durdaysL = (posixtime(stimmod.StimEnd(i)) - posixtime(stimmod.StimStart(i)))/(60*60*24);
            lowtime = cat(1,lowtime, durdaysL);
            durdaysL = [];
        elseif stimmod.Frequency(i) >= 50
            durdaysH = (posixtime(stimmod.StimEnd(i)) - posixtime(stimmod.StimStart(i)))/(60*60*24);
            hightime = cat(1,hightime, durdaysH);
            durdaysH = [];
        end 
    end 
    % calculate totals
    totaldayslowstim(s) = sum(lowtime);
    totaldayshighstim(s) = sum(hightime);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % find EMA scores during high and low freq stim
    lowstimEMA = [];
    highstimEMA = [];
    for i = 1: length(deltaselection.time)
        
        for ii = 1: height(stim)
            if deltaselection.time(i) > posixtime(stim.StimStart(ii)) && deltaselection.time(i) < posixtime(stim.StimEnd(ii))
                if stim.Frequency(ii) < 50
                    lowstimEMA = cat(1,lowstimEMA, deltaselection.signal(i));
                elseif stim.Frequency(ii) >= 50 
                    highstimEMA = cat(1,highstimEMA, deltaselection.signal(i));
                else
                    continue
                end
            else
                continue            
            end 
        end 
    end 
    totEMAscoreduringhighstim(s) = numel(highstimEMA);
    totEMAscoreduringlowstim(s) = numel(lowstimEMA);
    medEMAscoreduringhighstim(s) = median(highstimEMA);
    medEMAscoreduringlowstim(s) = median(lowstimEMA);
    meanEMAscoreduringhighstim(s) = mean(highstimEMA);
    meanEMAscoreduringlowstim(s) = mean(lowstimEMA);
    wilcoxonP_EMAscoreLowvHighstim(s) = ranksum(lowstimEMA, highstimEMA);  % caclulate wilcoxon rank sum

end 

%% Figure - EMA timeseries with recon, IED (spike) timeseries and seizures, IED and mood power spectra with noise floor
subjects = {'S1'};
Ns = [2000];
percent2drop = 10;
numrepeats = 10;

for i = 1: length(subjects)
    
    subject = subjects{i};
    N = Ns(i);
    
    % Load reference data from parameter sweep (N, dt, data, etc.) + sighits
    load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));

    % Load spike data 
    load(strcat(dest,subject,'_spikeinfo.mat'))
     
    fsize = 15;
    
    fig_timeseries = figure('Position',[10 10 1500 500])
    suptitle(subject)
    % plot the spike timeseries AND seizures
    subplot(2,6,[7,8,9,10])
    plot(spikeinfo.spikedates, spikeinfo.signal, 'Color','blue','DisplayName','IED rate')
%     title('Spike data')
%     xlabel('Days')
    hold on 
    plot(spikeinfo.szdates, max(spikeinfo.signal).*ones(1,length(spikeinfo.sz)),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','sz') % orange
    ylabel('IED rate')
    xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
    ylim([0, 1.1*max(spikeinfo.signal)])
    xl = xticklabels;
    for t= 1: length(xl)
        xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
    end 
    set(gca, 'xticklabel',xlnew)
    legend('Location','southeast')
    set(gca,'FontSize',fsize)
    
    % plot the mood signal reconstruction
    subplot(2,6,[1,2,3,4])
    plot(deltaselection.recondates, deltaselection.fig.reconsig,'Color','black','DisplayName','Recon')
    hold on
    plot(deltaselection.moodsampledates,deltaselection.signal, '*','Color','magenta','DisplayName','IMS')
    ylabel('EMA score')
    xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
    ylim([-36,36])
%     ylim([min(deltaselection.signal)-5, max(deltaselection.signal)+5])
    xl = xticklabels;
    for t= 1: length(xl)
        xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
    end 
    set(gca, 'xticklabel',xlnew)
    set(gca,'FontSize',fsize)
    if strcmp(subject,'M1')==1
        legend('Location','northeast')
    else
        legend('Location','southeast')
    end
    
    subplot(2,6,[5,6,11,12])    % overlay the BPDN and spike spectra 
    % spectrum for BPDN
    yyaxis right
    % no rescaling
    for iii = 1:size(deltaselection.fig.xreshuffled,1)
        plot(deltaselection.desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[1 0 1 0.1])
        hold on
    end 
    plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2,'-', 'Color','black', 'LineWidth',2);
    hold on
    plot(deltaselection.desiredperioddays(deltaselection.fig.sighits),deltaselection.fig.x(deltaselection.fig.sighits).^2,'o','Color','black','MarkerSize',7)
    % xlabel('Period (days)')
%     ylim([-100, 3*10^5])
    ylabel('Power EMA timeseries (BPWP)')
    set(gca,'FontSize',fsize)
    ax = gca;
    ax.YAxis(2).Color = 'black';
    % spectrum for cwt
    yyaxis left
    plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg, 'Color','blue', 'LineWidth',2);
    xlabel('Period (days)')
    ylabel({'';'';'Average power IED rate (CWT)'})
%     title('BPDN and CWT spectra')
    xlim([-2 175])
%     ylim([-100, 3*10^5])
    set(gca,'FontSize',fsize)
    ax = gca;
    ax.YAxis(1).Color = 'blue';
end

%% Figure 1: SCHEMATIC, FOR PAPER

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of the concurrent behavioral-electrophysiological timeseries
subject = 'S1';
percent2drop = 10;
numrepeats = 10;
N = 2000;

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));

% Load spike data 
load(strcat(dest, subject,'_spikeinfo.mat'))

fsize = 10;

fig_timeseries = figure('Position',[10 10 900 300])

% plot the mood signal samples
subplot(3,4,[1,2,3,4])
plot(deltaselection.moodsampledates,deltaselection.signal, '.','Color','magenta','MarkerSize',8,'DisplayName','IMS')
xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
ylim([-36,36])
set(gca,'FontSize',fsize)
legend('Location','northwest')
set(gca,'xticklabel',{[]})

% Plot the seizures when they occurred
subplot(3,4,[5,6,7,8])
plot(spikeinfo.szdates, ones(1,length(spikeinfo.sz)),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','sz') % orange
%ylabel('IED rate')
xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
ylim([0 2])
legend('Location','northwest')
set(gca,'FontSize',fsize)
set(gca,'xticklabel',{[]})

% plot the spike timeseries AND seizures
subplot(3,4,[9,10,11,12])
plot(spikeinfo.spikedates, spikeinfo.signal, 'Color','blue','DisplayName','IED rate')
xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
ylim([0, max(spikeinfo.signal)+20])
% remove the dates from the xticklabels - HIPAA
xl = xticklabels;
for t= 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
legend('Location','northwest')
set(gca,'FontSize',fsize)

%% FIGURE 2: FOR PAPER. Schematic figure 1/2 diagramming out BPWP outputs

subject = 'S1';
N = 2000; % optimal N for that subject
percent2drop = 10;
numrepeats = 10;

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

% reconstruct and pull out polynomial trend
maxdegree = deltaselection.maxdegree;
scalefactor = deltaselection.scalefactor;
f_k = deltaselection.f_k;
x = deltaselection.fig.x;
z = deltaselection.fig.z;
t_full = deltaselection.t_full;
recon_osc = idct_custom(f_k, x, scalefactor); % custom function to do idct based on f_k, matlab idct would give different output because it uses default basis
polycomp = [];
for i = 1:(maxdegree+1)
    polycomp(i,:) = t_full.^(i-1);
end 
polytrend = sum(z.*polycomp);
reconsig = recon_osc + sum(z.*polycomp);

% MAKE FIGURE -----------------------
fsize = 10;
fig_timeseries = figure('Position',[10 10 700 500])

% plot the mood signal samples
subplot(3,6,[1,2,3,4,5,6])
plot(deltaselection.moodsampledates,deltaselection.signal, '.','Color','magenta','MarkerSize',8,'DisplayName','EMA')
xlim([min(deltaselection.recondates) max(deltaselection.recondates)])
ylim([-40 35])
% remove the dates from the xticklabels - HIPAA
xl = xticklabels;
for t= 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')
% set(gca,'xticklabel',{[]})

% plot the spectrum (without significant hits)
subplot(3,6,[7,8])
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1.5);
set(gca,'FontSize',fsize)
xlim([-5,160])
ylim([-5000,3*10^5])

% plot the oscillations
subplot(3,6,[9,10])
plot(t_full, recon_osc, 'Color','black', 'LineWidth',0.1);
ylim([-50, 50])
set(gca,'xticklabel',{[]})

% plot the polynomial trend
subplot(3,6,[11,12])
plot(t_full, polytrend, 'Color','black', 'LineWidth',1.5);
set(gca,'xticklabel',{[]})

% plot the reconstruction
subplot(3,6,[13,14,15,16,17,18])
plot(deltaselection.recondates,deltaselection.fig.reconsig,'Color','black','DisplayName','Recon');
hold on
plot(deltaselection.moodsampledates,deltaselection.signal, '.','Color','magenta','MarkerSize',8,'DisplayName','EMA')
xlim([min(deltaselection.recondates) max(deltaselection.recondates)])
ylim([-55, 50])
% remove the dates from the xticklabels - HIPAA
xl = xticklabels;
for t= 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')
% set(gca,'xticklabel',{[]})

% plot the zoomed in reconstruction
fig = figure('Position',[10 10 300 150]);
plot(deltaselection.recondates,deltaselection.fig.reconsig,'Color','black','DisplayName','Recon');
hold on
plot(deltaselection.moodsampledates,deltaselection.signal, '*','Color','magenta','DisplayName','EMA')
xlim([deltaselection.recondates(1100) deltaselection.recondates(1350)])
ylim([-40, 30])
% remove the dates from the xticklabels - HIPAA
xl = xticklabels;
for t= 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)

%%%%%%%%%%%%%%%%%%%%%%%%
% subfigure explaining the significance evaluation/ significant peaks
fig_timeseries = figure('Position',[10 10 175 500])

% plot spectrum
subplot(3,1,1)
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1.5);
set(gca,'FontSize',fsize)
xlim([-5,160])
ylim([-5000,4*10^5])

% plot the spectrum noise floor
subplot(3,1,2)
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1.5);
for iii = 1:size(deltaselection.fig.xreshuffled,1)
    plot(deltaselection.desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.5 0.5 0.5 0.1])
    hold on
end 
set(gca,'FontSize',fsize)
xlim([-5,160])
ylim([-5000,4*10^5])

% plot sig stars and spectrum and noise floor 
subplot(3,1,3)
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1.5);
for iii = 1:size(deltaselection.fig.xreshuffled,1)
    plot(deltaselection.desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.5 0.5 0.5 0.1])
    hold on
end 
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1.5);
hold on
plot(deltaselection.desiredperioddays(sighits),deltaselection.fig.x(sighits).^2,'o','Color','black', 'MarkerSize',5)%,'LineWidth',1)
set(gca,'FontSize',fsize)
xlim([-5,160])
ylim([-5000,4*10^5])

%% FIGURE 3: FOR PAPER- Spike cycle methodology

subject = 'S1';
N = 2000; % optimal N for subject
percent2drop = 10;
numrepeats = 10;

%%%%%%%%%%%%%%%%%%%%%%
% figure of spike timeseries and cwt spectrum
    
% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

fsize = 10;
    
fig_timeseries = figure('Position',[10 10 800 200])
% plot the spike timeseries AND seizures
subplot(1,6,[1,2,3,4])
plot(spikeinfo.spikedates, spikeinfo.signal, 'Color','blue','DisplayName','IED rate')
hold on
plot(spikeinfo.szdates, ones(1,length(spikeinfo.szdates)).*max(spikeinfo.signal),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','seizure') % orange
xlim([min(spikeinfo.spikedates) max(spikeinfo.spikedates)])
ylim([0 max(spikeinfo.signal)*1.1])
xl = xticklabels;
for t= 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')

% plot spike spectrum
subplot(1,6,[5,6])
plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg, 'Color','blue', 'LineWidth',1);
set(gca,'FontSize',fsize)

% spike spectrum daily end zoom-in
% plot spike spectrum
fig_1 = figure('Position',[10 10 125 125])
plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg, 'Color','blue', 'LineWidth',1);
set(gca,'FontSize',fsize)

% filter hilbert output example
cycleranges ={[35,25]};% {[35,25]};%{[1.005,0.995];[10,5];[25,15];[35,25]}; % in days
perioddaycutoff = 2;

% Get phase at diff ranges per Leguia 
fs = 1 / (spikeinfo.time(10) - spikeinfo.time(9));
phase = [];
filtereds = [];
for f  = 1: length(cycleranges)

    % filter signal in range of interest
    [b,a] = butter(2,[(1./(cycleranges{f}.*24*60*60)./(fs/2))],'bandpass');
    y = filtfilt(b,a,spikeinfo.signal);
    filtereds(f,:) = y;
    % run hilbert 
    h = hilbert(y);
    phase(f,:) = angle(h); 
end 

% get indices of seizure events
szind = [];
for i = 1: length(sz)
    [val,loc] = min(abs(spikeinfo.time - sz(i)));
    szind(i) = loc;
end 

% plot timeseries of filtered signal (overlayed on raw) and phase of signal
fig = figure('Position',[10 10 800 200])
subplot(2,6,[1,2,3,4])
plot(spikeinfo.spikedates, spikeinfo.signal, 'Color','blue','DisplayName','IED rate')
% hold on
% plot(szdates, ones(1,length(spikeinfo.sz)).*max(spikeinfo.signal),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','sz') % orange
hold on
plot(spikeinfo.spikedates,real(h)+median(spikeinfo.signal),'Color','cyan','LineWidth',2,'DisplayName','filtered')
xlim([min(spikeinfo.spikedates) max(spikeinfo.spikedates)])
ylim([0 max(spikeinfo.signal)*1.1])
xl = xticklabels;
for t = 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')

% plot phase and seizures on top
subplot(2,6,[7,8,9,10])
plot(spikeinfo.spikedates, phase,'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName','phase') % color is dark purple
hold on
plot(spikeinfo.spikedates(szind), phase(szind),'o','LineWidth',1,'DisplayName','seizure')
xlim([min(spikeinfo.spikedates) max(spikeinfo.spikedates)])
ylim([-4 4])
xl = xticklabels;
for t = 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')

% Plot seizure phases in polar plot
fsize = 20;
fig = figure('Position',[10 10 400 400])
nbins = 24;%36;
polarhistogram(phase(:,szind),nbins,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.7)
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.ThetaDir = 'clockwise';
pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
pax.ThetaZeroLocation = 'left';
set(gca,'FontSize',fsize)

% Plot average seizure phase as mean vector 
[mean_theta ul ll] = circ_mean(phase(:,szind), [],2); % TRY wiht and without nbins as input
mean_rho = circ_r(phase(:,szind), [],[],2);
fig = figure('Position',[10 10 400 400])
polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',[0.8500 0.3250 0.0980])
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.ThetaDir = 'clockwise';
pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
pax.ThetaZeroLocation = 'left';
set(gca,'FontSize',fsize)

%% Figure 4: IED timeseries and seizure cycles per participant

% Figure spike timeseries per patient
% Figure classic seizure cycles per subject - common central frequencies -> physician confirmed seizures 
subjects = {'S1'};
Ns = [2000]; % optimal N per subject 
percent2drop = 10;
numrepeats = 10;
fsize=20;
    
for i = 1: length(subjects)
    
    subject = subjects{i};
    N=Ns(i);
    
    % Load reference data from parameter sweep (N, dt, data, etc.) + sighits
    load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
    % Load spike data
    load(strcat(dest, subject,'_spikeinfo.mat'))

    fs = 1 / (spikeinfo.time(2)-spikeinfo.time(1));
    
    % individual seizure phase plots
    polarcolor = [0.8500 0.3250 0.0980]; % orange'black';%[0 0.4470 0.7410].*0.5;
    alph = 0.8;
    nbins = 24;
    
    % intro processing
    cycleranges = {[0.8,1.2];[5,10];[15,25];[25,35]}; % in days - common cycles in focal epilepsies - can change to whichever cycles of interest
    [filtereds, phase] = firls_nick(spikeinfo.signal, fs, cycleranges);
    [szphases] = MoodCycles_getEventPhase_Nobuiltinfilter(spikeinfo.sz,spikeinfo.time,phase);

    %%%%%%%%%%%%%%%%%%%%%%
    % figure of spike timeseries and cwt spectrum

    fsize = 10;

    fig_timeseries = figure('Position',[10 10 1200 250])
    % plot the spike timeseries AND seizures
    subplot(1,6,[1,2,3,4,5])
    plot(spikeinfo.spikedates, spikeinfo.signal, 'Color',[0 0 1 0.1],'linewidth',1,'DisplayName','IED rate') %'linewidth',0.0001
    hold on
    plot(spikeinfo.szdates, ones(1,length(spikeinfo.sz)).*max(spikeinfo.signal),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','seizure') % orange
    xlim([min(spikeinfo.spikedates) max(spikeinfo.spikedates)])
    ylim([0 max(spikeinfo.signal)*1.1])
    xl = xticklabels;
    for t= 1: length(xl)
        xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
    end 
    set(gca, 'xticklabel',xlnew)
    set(gca,'FontSize',fsize)
    set(gca,'FontSize',fsize)

    % plot spike spectrum
    subplot(1,6,[6])
    plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg, 'Color','blue', 'LineWidth',1);
    xlim([-5 140])
    %ylim([-100, 3*10^5])
    set(gca,'FontSize',fsize)

    % POLAR PLOTS %%%%%%%%%%%%%%%%%
    fsize=20;
    fig = figure('Position',[10 10 2000 500])
    hours = hour(spikeinfo.szdates);
    % convert hours to phase
    hoursdeg = hours.*15; % convert to degrees 
    % one day phase plot
    subplot(1,4,1)
    polarhistogram(deg2rad(hoursdeg), deg2rad(0:15:360),'FaceColor',polarcolor,'FaceAlpha',alph)
    pax = gca;
    pax.ThetaTick = [90 270];
    pax.ThetaTickLabel = {'12AM';'12PM'};
    pax.ThetaZeroLocation = 'right'; % from playing around, very top is midnight, bottom is noon. 
    alpha = deg2rad(hoursdeg); % convert to radians
    [pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
    %   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
    if pval < 0.001
        titlestring = strcat('1 day',{' '}, '***');
    elseif pval < 0.01
        titlestring = {strcat('1 day cycle');strcat('Omnibus test, p = **')};
    elseif pval < 0.05
        titlestring = {strcat('1 day cycle');strcat('Omnibus test, p = *')};
    elseif pval >= 0.05
        titlestring = {strcat('1 day cycle');strcat('Omnibus test, p = ns')};
    end 
    title(titlestring)
    set(gca,'FontSize',fsize-2)

    % one week phase plot
    subplot(1,4,2)
    f=2;
    polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',alph)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    alpha = (szphases(f,:).*pi) ./180; % convert to radians
    [pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
    %   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
    if pval < 0.001
        titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
    elseif pval < 0.01
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = **')};
    elseif pval < 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = *')};
    elseif pval >= 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = ns')};
    end 
    title(titlestring)
    thetaticks(0:pi/2:2*pi)
    set(gca,'FontSize',fsize-2)

    % two-three week phase plot
    subplot(1,4,3)
    f=3;
    polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',alph)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    alpha = (szphases(f,:).*pi) ./180; % convert to radians
    [pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
    %   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
    if pval < 0.001
        titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
    elseif pval < 0.01
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = **')};
    elseif pval < 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = *')};
    elseif pval >= 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = ns')};
    end 
    title(titlestring)
    thetaticks(0:pi/2:2*pi)
    set(gca,'FontSize',fsize-2)

    % four week phase plot
    subplot(1,4,4)
    f=4;
    polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',alph)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    alpha = (szphases(f,:).*pi) ./180; % convert to radians
    [pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
    %   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
    if pval < 0.001
        titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
    elseif pval < 0.01
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = **')};
    elseif pval < 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = *')};
    elseif pval >= 0.05
        titlestring = {strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = ns')};
    end 
    title(titlestring)
    thetaticks(0:pi/2:2*pi)
    set(gca,'FontSize',fsize-2)

end

%% FIGURE 5: FOR PAPER - just showcasing the cycle reconstruction for each participant
%Note: all EMA samples are the same color
%Shows: IMS BPWP spectrum plust reconstructed timeseries with original timepoints

subject = 'S1';
N = 2000; % m1: 2000, m4:1500, m5:1500
side = 'left';%'right';
percent2drop = 10;
fsize=12;

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

% PLOT MOOD TIMESERIES recon with IMS samples 
fig = figure('Position',[10 10 800 200])
plot(deltaselection.recondates, deltaselection.fig.reconsig,'Color',[.7 .7 .7],'DisplayName','Recon')
hold on
plot(deltaselection.moodsampledates,deltaselection.signal, '*','Color','magenta','DisplayName','EMA')
ylabel('EMA score')
xlabel('Months')
xlim([recondates(1),recondates(end)])
ylim([-50, 50])
set(gca,'FontSize',fsize)
legend('Location','northeastoutside')
xl = xticklabels;
for t = 1: length(xl)
    xlnew{t} = cell2mat(regexp(xl{t},'[a-z]\w+','match','ignorecase'));
end 
set(gca, 'xticklabel',xlnew)
xlnew = [];
set(gca,'FontSize',fsize)

% PLOT BPWP SPECTRUM FOR SUBJECT - no thresholding, just sigstars 
fig2 = figure('Position',[10 10 300 250])
for iii = 1:size(deltaselection.fig.xreshuffled,1) %plot noisefloor 
    plot(deltaselection.desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.5 0.5 0.5 0.1])
    hold on
end 
plot(deltaselection.desiredperioddays(sighits),deltaselection.fig.x(sighits).^2,'o','Color','black','MarkerSize',7)
hold on
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2,'-', 'Color','black', 'LineWidth',2);
xlim([-5 170])
ylim([0 2*10^5]) % representing 5*10^4 for M1, 4*10^5 for rest
ylabel('Power')
xlabel('Period (days)')
set(gca,'FontSize',fsize)
ax = gca;

%% Figure 6. high+low+seizure+recon timeseries, band selected spectra, associated polar plots (mean vectors)

subject = 'S1';
N = 2000; % optimal N for that subject
percent2drop = 10;
numrepeats = 10;
percentilethreshold = 99;%
threshold = 3;% for high vs low events (tertiles)

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

[lowind, highind] = MoodCycles_getHighLowEvents(threshold,deltaselection.signal,deltaselection.time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE RANGES FOR BANDPASS BASED ON THRESHOLDED SPECTRAL CONCURRENCE *

p = prctile((deltaselection.fig.x.^2),percentilethreshold); % find threshold that is above 99th percentile of peaks overall
ab = find(deltaselection.fig.x.^2' > p);% find BPWP periods that are above p percentile, % ALSO have a minimum threshold in case sig peaks above percentile are still in noisefloor
besti = intersect(ab,deltaselection.fig.sighits);% then confirm that suprathresholds are also in sighits, only work with those moving forward
% best = indices of the bands to work with. 
cp = deltaselection.desiredperioddays(besti);% periods of the bands to work with - central frequencies basically
toremove = find(abs(diff(cp)) < 1); % if two central frequencies/periods are within one day of each other, omitone of them
cp(toremove) = [];
cp = fliplr(cp);
lp = cp - 0.1*cp;
up = cp + 0.1*cp;% use ps to define bandpass windows
% subsequent loops (filtering, figures, etc. will be based on length of cp)

numband = length(cp);
numcol = numband*2;
numrow = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsize=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MOOD TIMESERIES WITH HIGH/LOW EVENTS + Seizures
fig = figure('Position',[10 10 750 150])
plot(deltaselection.recondates, deltaselection.fig.reconsig,'Color',[0.7 0.7 0.7],'DisplayName','Recon') % color is gray
hold on
ylabel('EMA score')
hold on% plotting just high low mood events
plot(deltaselection.moodsampledates(highind),deltaselection.signal(highind),'*','Color',[0.4660 0.6740 0.1880],'DisplayName','High IMS sample') %green dark
hold on
plot(deltaselection.moodsampledates(lowind),deltaselection.signal(lowind),'*','Color',[0 0.4470 0.7410],'DisplayName','Low IMS sample') %dark blue
hold on 
plot(spikeinfo.szdates, 46.*ones(1,length(spikeinfo.szdates)),'o','Color',[0.8500 0.3250 0.0980],'DisplayName','Seizure') % orange
xlim([deltaselection.recondates(1),deltaselection.recondates(end)])
ylim([-50, 50])
legend('Location','northeastoutside')
xl = xticklabels;
for t = 1: length(xl)
    xlnew(t) = regexp(xl{t},'[a-z]\w+','match','ignorecase');
end 
set(gca, 'xticklabel',xlnew)
set(gca,'FontSize',fsize-2)
ylabel('EMA score')
xlabel('Months')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT BPWP SPECTRUM FOR SUBJECT - with thresholding - bars over highest
% amp trends + noise floor
fig = figure('Position',[10 10 250 200])
% plot noisefloor
for iii = 1:size(deltaselection.fig.xreshuffled,1) %plot noisefloor 
    plot(deltaselection.desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.5 0.5 0.5 0.1])
    hold on
end 
% spectrum for BPDN
ymax = 3*max(max(deltaselection.fig.xreshuffled.^2));
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2,'-', 'Color','black', 'LineWidth',2);
hold on
plot(deltaselection.desiredperioddays(deltaselection.fig.sighits),deltaselection.fig.x(deltaselection.fig.sighits).^2,'o','Color','black','MarkerSize',7)
% xlim([-10 170])
% ylim([0 2*10^5])
hold on
% plot shades for bandpass filters
for ii = 1:length(cp)
    l1 = [lp(ii) up(ii) up(ii) lp(ii)];
    l2 = [0 0 ymax ymax];
    patch(l1,l2,'k','FaceAlpha',0.2,'EdgeColor','none')
    hold on
end
ylabel('Power')
xlabel('Period (days)')
set(gca,'FontSize',fsize)
ax = gca;

% PLOT BPWP SPECTRUM FOR SUBJECT - with thresholding - bars over highest
% amp trends + IED avg power CWT
figx = figure('Position',[10 10 350 200])
% spectrum for BPWP
yyaxis right
ymax = 1.1*max(deltaselection.fig.x.^2);
ymax=2e5;
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1);
hold on
plot(deltaselection.desiredperioddays(deltaselection.fig.sighits),deltaselection.fig.x(deltaselection.fig.sighits).^2,'o','Color','black')
hold on
%yline(p,'--') % plots threshold used
% plot shades for bandpass filters
for ii = 1:length(cp)
    l1 = [lp(ii) up(ii) up(ii) lp(ii)];
    l2 = [0 0 ymax ymax];
    patch(l1,l2,'k','FaceAlpha',0.2,'EdgeColor','none')
    hold on
end
ylim([0,ymax])
xlim([-10 200])
ylabel('Power EMA')
set(gca,'FontSize',fsize)
ax = gca;
ax.YAxis(2).Color = 'black';
% spectrum for cwt
yyaxis left
plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg,'Color','blue','LineWidth',2)
xlabel('Period (days)')
ylabel('Power IED rate')
set(gca,'FontSize',fsize)
ax = gca;
ax.YAxis(1).Color = 'blue';

% PLOT BPWP SPECTRUM FOR SUBJECT - with thresholding - bars over highest
% amp trends + IED avg power CWT ZOOM IN HIGH FREQ END
figy = figure('Position',[10 10 400 200])
% spectrum for cwt
yyaxis left
plot(spikeinfo.perioddays_cwt, spikeinfo.meanpowerCWTreg,'Color','blue','LineWidth',2)
xlabel('Period (days)')
ylabel('Power IED rate (CWT)')
xlim([0 10])
set(gca,'FontSize',fsize)
ax = gca;
ax.YAxis(1).Color = 'blue';
% spectrum for BPWP
yyaxis right
ymax = 1.1*max(deltaselection.fig.x.^2);
plot(deltaselection.desiredperioddays, deltaselection.fig.x.^2, 'Color','black', 'LineWidth',1);
hold on
plot(deltaselection.desiredperioddays(deltaselection.fig.sighits),deltaselection.fig.x(deltaselection.fig.sighits).^2,'o','Color','black')
hold on
%yline(p,'--') % plots threshold used
% plot shades for bandpass filters
for ii = 1:length(cp)
    l1 = [lp(ii) up(ii) up(ii) lp(ii)];
    l2 = [0 0 ymax ymax];
    patch(l1,l2,'k','FaceAlpha',0.2,'EdgeColor','none')
    hold on
end
ylim([0,ymax])
ylabel('Power EMA (BPWP)')
set(gca,'FontSize',fsize)
ax = gca;
ax.YAxis(2).Color = 'black';

% PLOT POLAR HISTOGRAMS PER SPIKE PARADIGM 
% For later filtering and event phase identification

temp = diff(spikeinfo.time);
fs = 1 /temp(10);
fsize=40;
% suptitle(strcat(subject,{' '},side))
% Loop through cps'
numrow = 6;
nbins = 36;
threshold=3;
for i = 1: length(cp)
    
    % Filter and get relevant phases - depends on filter. for super slow
    % periods > 40 days, use filip
    if cp(i) > 40
        [filtered] = MoodSpikes_FilterFilip(spikeinfo.signal, spikeinfo.time,[up(i) lp(i)], []);
        phase = angle(hilbert(filtered));
    else 
        [filtered, phase] = firls_nick(spikeinfo.signal, fs, {[lp(i) up(i)]});
    end
    
    % Seizure events and self-events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_sz  = figure('Position',[10 10 800 800]) % 1800 1400
    [szphases] = MoodCycles_getEventPhase_Nobuiltinfilter(spikeinfo.sz,spikeinfo.time,phase);
    polarhistogram(szphases,nbins,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.8)
    hold on
%     polarhistogram(selfphases,nbins,'FaceColor','green','FaceAlpha',0.5)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaZeroLocation = 'left';
    pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
    pax.ThetaDir = 'clockwise';
    pax.RAxisLocation = pi;
    title(strcat(subject,{' '},num2str(round(lp(i),1)),'-',num2str(round(up(i),1)),{' '},'d IED cycle'))
    set(gca,'FontSize',fsize)
    
    % High/low mood resultant vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_highlow  = figure('Position',[10 10 800 800]) % 1800 1400
    [highphases, lowphases, allval_phases, random_highphases, random_lowphases, lowsamps, lowtimes, highsamps, hightimes] = MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter(phase, spikeinfo.time, deltaselection.signal, deltaselection.time, threshold);
    
    % STATS COMPARING MEAN PHASES
    [pval_lowVall, table] = circ_wwtest(lowphases, allval_phases);%, [w1, w2])
    [pval_highVall, table] = circ_wwtest(highphases, allval_phases);%, [w1, w2])
    [pval_lowVhigh, table] = circ_wwtest(highphases, lowphases);%, [w1, w2])
    if pval_lowVall < 0.001
        LvA = '***';
    elseif pval_lowVall < 0.01
        LvA = '**';
    elseif pval_lowVall < 0.05
        LvA = '*';
    elseif pval_lowVall >= 0.05
        LvA = 'ns';
    end
    if pval_highVall < 0.001
        HvA = '***';
    elseif pval_highVall < 0.01
        HvA = '**';
    elseif pval_highVall < 0.05
        HvA = '*';
    elseif pval_highVall >= 0.05
        HvA = 'ns';
    end
    if pval_lowVhigh < 0.001
        LvH = '***';
    elseif pval_lowVhigh < 0.01
        LvH = '**';
    elseif pval_lowVhigh < 0.05
        LvH = '*';
    elseif pval_lowVhigh >= 0.05
        LvH = 'ns';
    end 
    
    statstring = strcat('LvA=',LvA,', HvA=',HvA,', LvH=',LvH);
    statstring2 = strcat('High v Low',{' '},LvH,', p=',num2str(round(pval_lowVhigh,4)));
    statstring_simple = LvH;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %LOWs
    [mean_theta ul ll] = circ_mean(lowphases, [],2); % TRY wiht and without nbins as input
    mean_rho = circ_r(lowphases, [],[],2);
%     % plot CI around the mean
%     polarplot([ul,ul],[0 mean_rho],'LineWidth',2,'Color','blue')
%     hold on
%     polarplot([ll,ll],[0 mean_rho],'LineWidth',2,'Color','blue')
%     hold on
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color','blue')
    hold on
    % HIGHS
    [mean_theta ul ll] = circ_mean(highphases, [],2); % TRY wiht and without nbins as input
%     % plot CI around mean
%     polarplot([ul,ul],[0 mean_rho],'LineWidth',2,'Color','red')
%     hold on
%     polarplot([ll,ll],[0 mean_rho],'LineWidth',2,'Color','red')
    hold on
    mean_rho = circ_r(highphases, [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color','red')
    hold on
%     % LOWS RANDOM
%     [mean_theta ul ll] = circ_mean(reshape(random_lowphases,[1,numel(random_lowphases)]),[],2);% plot outline of low noisefloor %
%     mean_rho = circ_r(reshape(random_lowphases,[1,numel(random_lowphases)]), [],[],2);
%     polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',[0 0 1 0.4])
%     hold on
%     % HIGHS RANDOM
%     [mean_theta ul ll] = circ_mean(reshape(random_highphases,[1,numel(random_highphases)]),[],2);% plot outline of low noisefloor %
%     mean_rho = circ_r(reshape(random_highphases,[1,numel(random_highphases)]), [],[],2);
%     polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',	[1 0 0 0.4])
    % ALL PHASES
    [mean_theta ul ll] = circ_mean(reshape(allval_phases,[1,numel(allval_phases)]),[],2);% plot outline of low noisefloor %
    mean_rho = circ_r(reshape(allval_phases,[1,numel(allval_phases)]), [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',	[0 0 0 0.4])
%     % SEIZURES 
%     [mean_theta ul ll] = circ_mean(szphases, [],2); % TRY wiht and without nbins as input
%     mean_rho = circ_r(szphases, [],[],2);
%     polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',[0.8500 0.3250 0.0980])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    rlim([0 0.6])
    title1 = 'Mean phase of';
    title2 = 'High + low EMA samples';
    title4 = strcat(subject,{' '},num2str(round(lp(i),1)),'-',num2str(round(up(i),1)),{' '},'d IED cycle');
    title3 = statstring2{:};
    title(title3)
    set(gca,'FontSize',fsize)

    % PLOT LOWS, HIGHS, AND MEAN SEIZURE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOWs
    [mean_theta ul ll] = circ_mean(lowphases, [],2); % TRY wiht and without nbins as input
    mean_rho = circ_r(lowphases, [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',5,'Color',[0 0.4470 0.7410])
    hold on
    % HIGHS
    [mean_theta ul ll] = circ_mean(highphases, [],2); % TRY wiht and without nbins as input
    hold on
    mean_rho = circ_r(highphases, [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',5,'Color',[0.4660 0.6740 0.1880])
    hold on
    % SEIZURES
    [mean_theta ul ll] = circ_mean(szphases, [],2); % TRY wiht and without nbins as input
    hold on
    mean_rho = circ_r(szphases, [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',5,'Color',[0.8500 0.3250 0.0980])% Orange 
    hold on
   
    % ALL PHASES
    [mean_theta ul ll] = circ_mean(reshape(allval_phases,[1,numel(allval_phases)]),[],2);% plot outline of low noisefloor %
    mean_rho = circ_r(reshape(allval_phases,[1,numel(allval_phases)]), [],[],2);
    polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',5,'Color',	[0 0 0 0.4])
%     % SEIZURES 
%     [mean_theta ul ll] = circ_mean(szphases, [],2); % TRY wiht and without nbins as input
%     mean_rho = circ_r(szphases, [],[],2);
%     polarplot([mean_theta, mean_theta],[0, mean_rho],'LineWidth',4,'Color',[0.8500 0.3250 0.0980])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    rlim([0 0.6])
    title1 = 'Mean phase of';
    title2 = 'High + low EMA samples';
    title4 = strcat(subject,{' '},num2str(round(lp(i),1)),'-',num2str(round(up(i),1)),{' '},'d IED cycle');
    title3 = statstring2{:};
    title(title3)
    set(gca,'FontSize',fsize)
        
    % All events phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_all  = figure('Position',[10 10 800 800]) % 1800 1400
    polarhistogram(allval_phases,nbins,'EdgeColor','black','DisplayStyle','stairs');% plot outline of low noisefloor %
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaZeroLocation = 'left';
    pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
    pax.ThetaDir = 'clockwise';

    alpha = (lowphases.*pi) ./180; % convert to radians
    [pval_k, k, K] = circ_kuipertest(lowphases, highphases);% res, vis_on);
    [pval_w, table] = circ_wwtest(lowphases, highphases); % 
    title1 = 'Phase dist. of';
    title2 = 'all EMA samples';
    title(title4{:})
    set(gca,'FontSize',fsize)
    
end 

%% FIGURE 7: FOR PAPER.  vectors overlayed on sz histogram: POlar histograms (non-averaged) from filtered cycles 

subject = 'S1';
N = 2000; % optimal N for subject 
percent2drop = 10;
numrepeats = 10;
percentilethreshold = 99;%

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

% For later filtering and event phase identification
temp = diff(spikeinfo.time);
fs = 1 /temp(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE RANGES FOR BANDPASS BASED ON THRESHOLDED SPECTRAL CONCURRENCE *

p = prctile((deltaselection.fig.x.^2),percentilethreshold); % find threshold that is above 99th percentile of peaks overall
ab = find(deltaselection.fig.x.^2' > p);% find BPWP periods that are above p percentile, % ALSO have a minimum threshold in case sig peaks above percentile are still in noisefloor
besti = intersect(ab,deltaselection.fig.sighits);% then confirm that suprathresholds are also in sighits, only work with those moving forward
% best = indices of the bands to work with. 
cp = deltaselection.desiredperioddays(besti);% periods of the bands to work with - central frequencies basically
toremove = find(abs(diff(cp)) < 1); % if two central frequencies/periods are within one day of each other, omitone of them
cp(toremove) = [];
cp = fliplr(cp);
lp = cp - 0.1*cp;
up = cp + 0.1*cp;% use ps to define bandpass windows
% subsequent loops (filtering, figures, etc. will be based on length of cp)

numband = length(cp);
numcol = numband*2;
numrow = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT POLAR HISTOGRAMS PER SPIKE PARADIGM
fsize=40;
% suptitle(strcat(subject,{' '},side))
% Loop through cps'
numrow = 6;
nbins = 36;
threshold=3;
for i = 1: length(cp)
    
    % Filter and get relevant phases - depends on filter. for super slow
    % periods > 40 days, use filip
    if cp(i) > 40
        [filtered] = MoodSpikes_FilterFilip(spikeinfo.signal, spikeinfo.time,[up(i) lp(i)], []);
        phase = angle(hilbert(filtered));
    else 
        [filtered, phase] = firls_nick(spikeinfo.signal, fs, {[lp(i) up(i)]});
    end
    
    % Seizure events 
    [szphases] = MoodCycles_getEventPhase_Nobuiltinfilter(spikeinfo.sz,spikeinfo.time,phase);
    [highphases, lowphases, allval_phases, random_highphases, random_lowphases, lowsamps, lowtimes, highsamps, hightimes] = MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter(phase, spikeinfo.time, deltaselection.signal, deltaselection.time, threshold);

    % STATS COMPARING MEAN PHASES
    [pval_lowVall, table] = circ_wwtest(lowphases, allval_phases);%, [w1, w2])
    [pval_highVall, table] = circ_wwtest(highphases, allval_phases);%, [w1, w2])
    [pval_lowVhigh, table] = circ_wwtest(highphases, lowphases);%, [w1, w2])
    if pval_lowVall < 0.001
        LvA = '***';
    elseif pval_lowVall < 0.01
        LvA = '**';
    elseif pval_lowVall < 0.05
        LvA = '*';
    elseif pval_lowVall >= 0.05
        LvA = 'ns';
    end
    if pval_highVall < 0.001
        HvA = '***';
    elseif pval_highVall < 0.01
        HvA = '**';
    elseif pval_highVall < 0.05
        HvA = '*';
    elseif pval_highVall >= 0.05
        HvA = 'ns';
    end
    if pval_lowVhigh < 0.001
        LvH = '***';
    elseif pval_lowVhigh < 0.01
        LvH = '**';
    elseif pval_lowVhigh < 0.05
        LvH = '*';
    elseif pval_lowVhigh >= 0.05
        LvH = 'ns';
    end 
    
    statstring = strcat('LvA=',LvA,', HvA=',HvA,', LvH=',LvH);
    statstring2 = strcat('High v Low',{' '},LvH,', p=',num2str(round(pval_lowVhigh,4)));
    statstring_simple = LvH;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EEG SZ
    % vectors overlayed on seizure histogram(to make polarplot and polar
    % histogram overlay work, reformat vector into a histogram bin of position theta and length rho)
    % NOTE vectors are rescaled so length is more aparent in setting of
    % histogram 
    %SZ
    linethickness = 6;%4;
    vectoralpha = 0.8;
    facealpha = 0.1;
    fig_szoverlay  = figure('Position',[10 10 800 800]) % 1800 1400
    h = polarhistogram(szphases,nbins,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',facealpha, 'EdgeColor',[0.8500 0.3250 0.0980],'EdgeAlpha',0.4);% plot outline of low noisefloor %
    rescalefactor = 2*max(h.Values);
    hold on
    % ALL PHASES
    [mean_theta ul ll] = circ_mean(reshape(allval_phases,[1,numel(allval_phases)]),[],2);% plot outline of low noisefloor %
    mean_rho = circ_r(reshape(allval_phases,[1,numel(allval_phases)]), [],[],2);
    polarhistogram('BinEdges',[mean_theta, mean_theta],'BinCounts',mean_rho*rescalefactor,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7],'EdgeAlpha',vectoralpha,'LineWidth',linethickness)
    %LOWs
    [mean_theta ul ll] = circ_mean(lowphases, [],2); % TRY wiht and without nbins as input
    mean_rho = circ_r(lowphases, [],[],2);
    polarhistogram('BinEdges',[mean_theta, mean_theta],'BinCounts',mean_rho*rescalefactor, 'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'EdgeAlpha',vectoralpha,'LineWidth',linethickness)
    hold on
    % HIGHS
    [mean_theta ul ll] = circ_mean(highphases, [],2); % TRY wiht and without nbins as input
    mean_rho = circ_r(highphases, [],[],2);
    polarhistogram('BinEdges',[mean_theta, mean_theta],'BinCounts',mean_rho*rescalefactor,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880],'EdgeAlpha',vectoralpha,'LineWidth',linethickness)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'left';
    set(gca,'FontSize',fsize)
    rticks(pax,[round(0.5*rescalefactor) rescalefactor])
    rticklabels(pax,{num2str(round(0.5*rescalefactor)),num2str(rescalefactor)})
    
end 

%% Tables for main results: table of circ stats results, stat outputs 

subject = 'S1';
N = 2000;% optimal N for subject
percent2drop = 10;
numrepeats = 10;
percentilethreshold = 99;%

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

% For later filtering and event phase identification
temp = diff(spikeinfo.time);
fs = 1 /temp(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE RANGES FOR BANDPASS BASED ON THRESHOLDED SPECTRAL CONCURRENCE *

p = prctile((deltaselection.fig.x.^2),percentilethreshold); % find threshold that is above 99th percentile of peaks overall
ab = find(deltaselection.fig.x.^2' > p);% find BPWP periods that are above p percentile, % ALSO have a minimum threshold in case sig peaks above percentile are still in noisefloor
besti = intersect(ab,deltaselection.fig.sighits);% then confirm that suprathresholds are also in sighits, only work with those moving forward
% best = indices of the bands to work with. 
cp = deltaselection.desiredperioddays(besti);% periods of the bands to work with - central frequencies basically
toremove = find(abs(diff(cp)) < 1); % if two central frequencies/periods are within one day of each other, omitone of them
cp(toremove) = [];
cp = fliplr(cp);
lp = cp - 0.1*cp;
up = cp + 0.1*cp;% use ps to define bandpass windows
% subsequent loops (filtering, figures, etc. will be based on length of cp)

numband = length(cp);
threshold=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats table 1
p_omn_sz = []; 
p_omn_high = [];
p_omn_low = []; % resultant
r_sz = [];
r_high = [];
r_low = [];

% stats table 2
pvals_lowVsz_2samp = [];
pvals_highVsz_2samp = [];
pvals_lowVhigh_2samp = [];
pvals_anova_ish = [];
warnings_ANOVA = {};
warnings_lowVsz = {};
warnings_highVsz = {};
warnings_lowVhigh = {};
pvals_watsonU2_lowVsz_2samp = [];
pvals_watsonU2_highVsz_2samp = [];
pvals_watsonU2_lowVhigh_2samp = [];
U2vals_watsonU2_lowVsz_2samp = [];
U2vals_watsonU2_highVsz_2samp = [];
U2vals_watsonU2_lowVhigh_2samp = [];
warnings_U2_lowVsz = {};
warnings_U2_highVsz = {};
warnings_U2_lowVhigh = {};

for i = 1: length(cp)
    
    % Filter and get relevant phases - depends on filter. for super slow
    % periods > 40 days, use filip
    if cp(i) > 40
        [filtered] = MoodSpikes_FilterFilip(spikeinfo.signal, spikeinfo.time,[up(i) lp(i)], []);
        phase = angle(hilbert(filtered));
    else 
        [filtered, phase] = firls_nick(spikeinfo.signal, fs, {[lp(i) up(i)]});
    end
    
    % Seizure events and self-events 
    [szphases] = MoodCycles_getEventPhase_Nobuiltinfilter(spikeinfo.sz,spikeinfo.time,phase);
    % High/low mood resultant vectors 
    [highphases, lowphases, allval_phases, random_highphases, random_lowphases, lowsamps, lowtimes, highsamps, hightimes] = MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter(phase, spikeinfo.time, deltaselection.signal, deltaselection.time, threshold);
    
    %%%%%%%%%%%%%%%%%%%%%%%% STATS TABLE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Omnibus test - Hodges-Ajne 
    % does sample distribution deviate from uniformity?
    p_omn_sz(i,1)= circ_otest(szphases);
    p_omn_high(i,1)= circ_otest(highphases);
    p_omn_low(i,1)= circ_otest(lowphases);
    
    % resultant vector lengths
    r_sz(i,1) = circ_r(szphases');
    r_high(i,1) = circ_r(highphases');
    r_low(i,1) = circ_r(lowphases');
    
    % %%%%%%%%%%%%%%%%%%%%%% STATS TABLE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Watson williams
    % circular analog of two sample t test, one way anova
    % assesses the question whether the mean directions of two or more
    % groups are identical or not? 
    
    % Setup data for ANOVA-style (are means diff?)
    alphas = [highphases, lowphases, szphases];
    idx = [1.*ones(length(highphases),1);2.*ones(length(lowphases),1);3.*ones(length(szphases),1)]';
    warning(''); % clears warning between iterations
    [pval, tab] = circ_wwtest(alphas, idx);
    pvals_anova_ish(i,1) = pval;
    warnings_ANOVA{i,1} = lastwarn;
    
    % Setup for two-sample style (are means diff?)
    warning(''); % clears warning between iterations
    [p_lowVsz] = circ_wwtest(lowphases, szphases);%, [w1, w2])
    pvals_lowVsz_2samp(i,1) = p_lowVsz;
    warnings_lowVsz{i,1} = lastwarn;
    
    warning(''); % clears warning between iterations
    [p_highVsz] = circ_wwtest(highphases, szphases);%, [w1, w2])
    pvals_highVsz_2samp(i,1) = p_highVsz;
    warnings_highVsz{i,1} = lastwarn;
    
    warning(''); % clears warning between iterations
    [p_lowVhigh] = circ_wwtest(highphases, lowphases);%, [w1, w2]
    pvals_lowVhigh_2samp(i,1) = p_lowVhigh;
    warnings_lowVhigh{i,1} = lastwarn;
    
    
    %%%%%%%%%%%%%%%%%%% Stats table 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Watsons U2 test 
    % % p: probability that the two samples come from the same distribution
    % U2_obs:   observed value of Watson's U2 statistic
    % from https://github.com/pierremegevand/watsons_u2 
    % best option for comparing means circular data: Advice on comparing two independent samples of circular data in biology
    % Nature: Landler et al. 2021
    % https://www.nature.com/articles/s41598-021-99299-5 
    
    warning(''); % clears warning between iterations
    [pU2_lowVsz,U2_lowSz] = watsons_U2_approx_p(lowphases', szphases');%, [w1, w2])
    pvals_watsonU2_lowVsz_2samp(i,1) = pU2_lowVsz;
    U2vals_watsonU2_lowVsz_2samp(i,1) = U2_lowSz;
    warnings_U2_lowVsz{i,1} = lastwarn;
    
    warning(''); % clears warning between iterations
    [pU2_highVsz, U2_highVsz] = watsons_U2_approx_p(highphases', szphases');%, [w1, w2])
    pvals_watsonU2_highVsz_2samp(i,1) = pU2_highVsz;
    U2vals_watsonU2_highVsz_2samp(i,1) = U2_highVsz; 
    warnings_U2_highVsz{i,1} = lastwarn;
    
    warning(''); % clears warning between iterations
    [pU2_lowVhigh,U2_lowVhigh] = watsons_U2_approx_p(highphases', lowphases');%, [w1, w2]
    pvals_watsonU2_lowVhigh_2samp(i,1) = pU2_lowVhigh;
    U2vals_watsonU2_lowVhigh_2samp(i,1) = U2_lowVhigh;
    warnings_U2_lowVhigh{i,1} = lastwarn;
        
end 

clear table

% stats table 1 (non/uniformity tests)
stats1 = table;
stats1.cycledays = cp';
stats1.p_omn_sz = p_omn_sz;
stats1.resultant_sz = r_sz;
stats1.p_omn_high = p_omn_high;
stats1.resultant_high = r_high;
stats1.p_omn_low = p_omn_low;
stats1.resultant_low = r_low;
writetable(stats1,strcat(dest,subject,'_omnibus_and_resultants_perperiod.csv'))

% stats table 1, Watson williams 
stats2 = table;
stats2.cycledays = cp';
stats2.p_anova = pvals_anova_ish;
stats2.warn_anova = warnings_ANOVA;
stats2.p_lowVsz = pvals_lowVsz_2samp;
stats2.warn_lowVsz = warnings_lowVsz;
stats2.p_highVsz = pvals_highVsz_2samp;
stats2.warn_highVsz = warnings_highVsz;
stats2.p_lowVhigh = pvals_lowVhigh_2samp;
stats2.warn_lowVhigh = warnings_lowVhigh;
writetable(stats2,strcat(dest,subject,'_watsonWilliamstwosamp_and_ANOVA_perperiod.csv'))

% stats table 2, watson U2 test
stats3 = table;
stats3.cycledays = cp';
stats3.p_lowVsz = pvals_watsonU2_lowVsz_2samp;
stats3.U2_lowVsz = U2vals_watsonU2_lowVsz_2samp;
stats3.warn_lowVsz = warnings_U2_lowVsz;
stats3.p_highVsz = pvals_watsonU2_highVsz_2samp;
stats3.U2_highVsz = U2vals_watsonU2_highVsz_2samp;
stats3.warn_highVsz = warnings_U2_highVsz;
stats3.p_lowVhigh = pvals_watsonU2_lowVhigh_2samp;
stats3.U2_lowVhigh = U2vals_watsonU2_lowVhigh_2samp;
stats3.warn_lowVhigh = warnings_U2_lowVhigh;
writetable(stats3,strcat(dest,subject,'_WatsonU2_perperiod.csv'))

%% SUPPLEMENTARY FIGURE 2-4: Mood vs stimulation parameters - same fig for frequency and amplitude
 
% did not simulate stim parameter table here, can input one's own
% necessary columns: subject, stimstart time, stimendtime, then columns for
% parameters (frequency, amplitude)

subjects = {'S1'};
Ns=[2000]; % optimal N for subject
percent2drop = 10;
numrepeats = 10;
    
for i = 1: length(subjects)
    
    subject = subjects{i};
    N=Ns(i);

    % Load reference data from parameter sweep (N, dt, data, etc.) + sighits
    load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
    % Load spike data
    load(strcat(dest, subject,'_spikeinfo.mat'))

    % Load stimulation parameters
    tbase = []; % first timestamp of perferred reference timeseries, usually first EMA sample
    stimfile = '';
    stim = readtable(stimfile);
    stimstarts = (stim.start_uutc - tbase) ./ (24*60*60); % converting to days
    stimends = (stim.stop_uutc - tbase) ./ (24*60*60);
   
    % Plot........ 
    fsize = 15;
    fig = figure('Position',[10,10,2000,1500])
    
    % STIM FREQUENCY
    subplot(3,4,[1,2,3,4])
    for s = 1:length(stimstarts)
        line([stimstarts(s) stimends(s)],[stim.Frequency(s) stim.Frequency(s)],'LineWidth',2,'Color','blue')
        hold on
    end 
    title({'';'';'';'Stimulation frequency'})
    xlabel('Days')
    ylabel('Frequency (Hz)')
    set(gca,'FontSize',fsize)
    xlim([0,max(deltaselection.moodrecondates)])
    ylim([-5,155])
    
    % plot the mood signal reconstruction
    subplot(3,4,[9,10,11,12])
    plot(deltaselection.moodrecondates, deltaselection.fig.reconsig,'Color','black')
    hold on
    plot(deltaselection.moodsampledates,deltaselection.signal, '*','Color','magenta')
    title('IMS scores and reconstructed signal')
    xlabel('Days')
    ylabel('IMS total score')
    xlim([0,max(deltaselection.moodrecondates)])
    set(gca,'FontSize',fsize)
    
    % STIM AMPLITUDE
    % plot the weekly seizures 
    subplot(3,4,[5,6,7,8])
    for s = 1:length(stimstarts)
        line([stimstarts(s) stimends(s)],[stim.Amps(s) stim.Amps(s)],'LineWidth',2,'Color','red')
        hold on
    end 
    title({'';'';'';'Stimulation amplitude'})
    xlabel('Days')
    ylabel('Amplitude (mA)')
    set(gca,'FontSize',fsize)
    xlim([0,max(deltaselection.moodrecondates)])
    ylim([-2, 10])
    
end

%% SUPPLEMENTARY FIGURE 5 - Data output: Power of IED cycle at BPWP dictated period of choice

subject = {'S1'};
N = [2000];
percent2drop = 10;
numrepeats = 10;
percentilethreshold = 99;%

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
load(strcat(dest,subject,'_Outputfor_','MSE_','percent2drop',num2str(percent2drop),'_N',num2str(N),'_traintest_forfig.mat'));
% Load spike data
load(strcat(dest, subject,'_spikeinfo.mat'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE RANGES FOR BANDPASS BASED ON THRESHOLDED SPECTRAL CONCURRENCE *

p = prctile((deltaselection.fig.x.^2),percentilethreshold); % find threshold that is above 99th percentile of peaks overall
ab = find(deltaselection.fig.x.^2' > p);% find BPWP periods that are above p percentile, % ALSO have a minimum threshold in case sig peaks above percentile are still in noisefloor
besti = intersect(ab,deltaselection.fig.sighits);% then confirm that suprathresholds are also in sighits, only work with those moving forward
% best = indices of the bands to work with. 
cp = deltaselection.desiredperioddays(besti);% periods of the bands to work with - central frequencies basically
toremove = find(abs(diff(cp)) < 1); % if two central frequencies/periods are within one day of each other, omitone of them
cp(toremove) = [];
cp = fliplr(cp);
lp = cp - 0.1*cp;
up = cp + 0.1*cp;% use ps to define bandpass windows
% subsequent loops (filtering, figures, etc. will be based on length of cp)

numband = length(cp);
numcol = numband*2;
numrow = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Position',[10 10 1200 1600]) 
percentiles = [];
% Histogram for each CP
for ii = 1:length(cp)

    subplot(1,length(cp),ii)

    % CWT Spec as histogram
    histogram(spikeinfo.meanpowerCWTreg,20,'FaceColor',[0.3010 0.7450 0.9330])
    [c ind] = min(abs(spikeinfo.perioddays_cwt-cp(ii)));
    poi = spikeinfo.meanpowerCWTreg(ind);
    label = strcat(num2str(round(cp(ii),2))," ",'d');
    xline(poi,":",label,'LabelOrientation','horizontal','LineWidth',2)
    ylabel('Count')
    xlabel('IED CWT spectral power')

    % calculate percentiles 
    centile = comp_percentile(spikeinfo.meanpowerCWTreg,poi);
    percentiles(ii,1) = centile;
end

% output percentiles table and save
clear table
percs = table;
percs.cycledays = cp';
percs.percentiles = percentiles;
savename = strcat(dest,subject,'_percentileOfBPWPCyleOfInterestWithinIEDCWTPowerDistrib.csv');
writetable(percs,savename{:})
