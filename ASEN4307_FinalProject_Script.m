%% ASEN4307 Final Project: GPS Signal, Ionospheric Impact, Solar Flux

%%
% *Lukas Zumwalt*
%
% *12/13/2024*

%% Context
%
% For this investigation, I am seeking to analyze and assess raw GPS data
% and parameters over time, specifically looking at ionospheric and
% tropospheric fluctuations and if they correlate with solar phenomena.
%
% For data, I am referencing publicly available NASA CDDIS broadcast
% ephemeris and SP3 acquisition data in order to understand and
% describe the observations and experience of GPS satellite signals.
%
% # https://cddis.nasa.gov/archive/gnss/products/
% 
% I am also using solar flux data available from the government of Canada
% to describe and understand any patterns present in that data set.
%
% # https://www.spaceweather.gc.ca/forecast-prevision/solar-solaire/solarflux/sx-5-en.php

%% Initialization and Data Loading
clc; clear; close all;

% Constants
c = 299792458;  % Speed of light (m/s)
omegaE = 7.2921151467e-5; % Earth's rotation rate (rad/s)
f1 = 1575.42e6; % GPS L1 frequency

% File paths
broadcastPath = 'data/';  % e.g., 'data/broadcast/'
precisePath = 'data/'; % e.g., 'data/precise/'
% rinexFile = 'data_rinex_NIST00USA_R_2024_239.mat';


% Dates for analysis
yyyy = 2004:2023;
year_sequence = linspace(2004,2023,80);


%% Solar Flux Data

% Define the input file name
inputFileName = 'data/flux_data.txt';

% Read the file as a table using fixed-width text format
opts = fixedWidthImportOptions('NumVariables', 7);

% Define variable names and types
opts.VariableNames = {'fluxdate', 'fluxtime', 'fluxjulian', 'fluxcarrington', ...
                      'fluxobsflux', 'fluxadjflux', 'fluxursi'};
opts.VariableTypes = {'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double'};

% Define the column widths based on the format in the file
opts.DataLines = [3 Inf]; % Skip the header lines (adjust as needed)
opts.VariableWidths = [10, 10, 15, 15, 13, 13, 11];

% Read the table
fluxTable = readtable(inputFileName, opts);
fluxTable = fluxTable(:,1:7);
% Prune NaN's (replace with 0)
fluxTable = fillmissing(fluxTable,'constant',0);


% Save the table into a MATLAB file
outputFileName = 'flux_data.mat';
save(outputFileName, 'fluxTable');

% Visualize raw data
figure;
grid on;
hold on;
plot(linspace(2004,2023,length(fluxTable.fluxobsflux)),fluxTable.fluxobsflux,'.-');
plot(linspace(2004,2023,length(fluxTable.fluxadjflux)),fluxTable.fluxadjflux,'.-');
plot(linspace(2004,2023,length(fluxTable.fluxursi)),fluxTable.fluxursi,'.-');
xlabel('Epoch');
ylabel('Flux');
title('Solar Flux over Time, Canada');
legend({'OBS','ADJ','URSI'});
axis tight

%%
%
% There are definite peaks in the available flux data, which completely
% dwarf the rest of the data in magnitude.  I am assuming these to be
% relative events of solar flares or hyperlocal chaotic periods in the
% sun's corona.

%% Data Inspection

% Solar Flux FFT
freq_fluxobs = abs(fft(fluxTable.fluxobsflux));

figure
hold on;
grid on;
plot(linspace(2004,2023,length(freq_fluxobs)),freq_fluxobs)
xlabel('Frequency')
ylabel('Amplitude')
title('Observed Flux FFT, Wide')

%%
%
% At this point, I can tell I need a low-pass filter, as the 0-frequency is
% much louder than the rest of the spectrum and makes it difficult to
% discern trends.  Since I am seeking more than the prevelance of
% 1-year-resolution trends, I can remove that and re-view.

% Solar Flux FFT, selecting specific band past 0-frequency to end of first
% half of frequencies.
freq_fluxobs = freq_fluxobs(10:floor(end/2));

figure
hold on;
grid on;
plot(abs(freq_fluxobs))
xlabel('Frequency')
ylabel('Amplitude')
title('0-Frequency Filtered, Observed Flux FFT')


%%
% Compare against a periodogram

figure
periodogram(fluxTable.fluxobsflux)


%% GPS Data
%
% In this circumstance, I am selecting a BRDC ephemeric report
% once every 30 days, continuing at a round modulo of 30. 
% (i.e. doy 30, 60, 90, 120, ... etc)

% Resolution arrays
yyyy = 2004:2023;
year_sequence = linspace(2004,2023,80);

% constant
c = 299792458;
omegaE = 7.2921151467e-5;
R_E = 6378.137e3;
f1 = 1575.42e6; f2 = 1227.6e6;

% NIST ECEF/LLA Location
selected_site_ecef = [-1288398.574, -4721696.936, 4078625.349];
selected_site_lla = ecef2lla(selected_site_ecef);

% Start timing data collection
% tic

% Ephemeris
select_doys = {'0900','1800','2700','3600'};
ephem_all = {0};
iono_all = {0};
brdc_iter = 1;
for year = yyyy
    brdc_yyyy = dir(['data/brdc/' num2str(year) '/*']);
    for ii = 1:length(brdc_yyyy)
        for jj = 1:length(select_doys)
            if contains(brdc_yyyy(ii).name,select_doys(jj))
                % By using |dir| attributes, we can avoid issues with different
                % file-name formats, which change over time.
                [ephem_all{brdc_iter},iono_all{brdc_iter}] = ...
                    read_GPSbroadcast([brdc_yyyy(ii).folder '\'  brdc_yyyy(ii).name]);
                % Increment
                brdc_iter = brdc_iter + 1;
            end
        end
    end
end

% SP3
sp3_all = {0};
sp3_iter = 1;
for year = yyyy
    sp3_yyyy = dir(['data/sp3/' num2str(year) '/*.sp3']);
    for ii = 1:length(sp3_yyyy)
        % By using |dir| attributes, we can avoid issues with different
        % file-name formats, which change over time.
        sp3_all{sp3_iter} = read_sp3([sp3_yyyy(ii).folder '\'  sp3_yyyy(ii).name]);
        sp3_iter = sp3_iter + 1;
    end
end

% Finish timing data collection
% toc

%%
% Orbits and Clock Terms

Az_all = {0};
El_all = {0};
Range_all = {0};
tow_sequence = {0};
gps_week_sequence = {0};
azelra_iter = 1;
constellation = 1; % 1 for GPS

% Iterate over each epoch, 4 quarters per year from 2004 to 2023
for epoch = 1:length(ephem_all)

    % Save to avoid function-return-value indexing
    sp3_epoch = sp3_all{epoch};

    % At this epoch, filter for each individual PRN
    for i = 1:32
        selected_prn = i;

        if size(sp3_epoch,2) == 8
            sp3_line_index = find((sp3_epoch(:,8)==constellation)&(sp3_epoch(:,3)==selected_prn));
        else
            sp3_line_index = find((sp3_epoch(:,3)==selected_prn));
        end

        % Grab some relative time arrays
        tow_sequence{epoch,selected_prn} = sp3_epoch(sp3_line_index,2);
        gps_week_sequence{epoch,selected_prn} = sp3_epoch(sp3_line_index,1);

        % Calculate ephemeris orbits & clock impacts
        [~,satPos_eph_all{epoch,selected_prn}, ...
            satVel_eph_all{epoch,selected_prn}, ...
            satClkCorr_eph_all{epoch,selected_prn}, ...
            relCorr_eph_all{epoch,selected_prn}, ...
            tgd_eph_all{epoch,selected_prn}] = ...
        eph2pvt(ephem_all{epoch},[gps_week_sequence{epoch,selected_prn},tow_sequence{epoch,selected_prn}],selected_prn);

        % Save position relative to NIST
        satPos_sp3_thisPRN = sp3_epoch(sp3_line_index,4:6)*1e3;
        [Az_all{epoch,selected_prn},El_all{epoch,selected_prn},Range_all{epoch,selected_prn}] = ...
            compute_azelrange(selected_site_ecef, satPos_sp3_thisPRN);
    end
end

%%
%
% At this point, I have three 80x32 (epochs x prn) cell arrays containing
% vectors for altitude, elevation, and range relative to the NIST
% observatory; as well as an 80x32 set of ephemeris-derived state
% characteristics of the satellites, such as position/velocity and clock
% and relativistic corrections.
%
% I also have a 1x80 cell array of both ephemeris and ionospheric
% parameters for each epoch, describing a daily figure of satellite orbital
% characteristics and atmospheric conditions.  I need this to generate a
% model of ionospheric delay experienced by the GPS signal travel.



%% Modeling Ionospheric Delay
%
% https://gssc.esa.int/navipedia/index.php/Klobuchar_Ionospheric_Model

I_s = {0};
for epoch = 1:length(ephem_all)
    for selected_prn = 1:32

        % Localize azimuth and elevation
        Az_selected_prn = Az_all{epoch,selected_prn};
        El_selected_prn = El_all{epoch,selected_prn};

        if isempty(Az_selected_prn) || isempty(El_selected_prn)
            I_s{epoch,selected_prn} = NaN;
            continue;
        end

        phi_E = 0.0137./(El_selected_prn/180+0.11)-0.022; % semicircle
        lat_I = selected_site_lla(1)/180+cosd(Az_selected_prn).*phi_E; % semicircle
        lat_I(lat_I>0.416) = 0.416;
        lat_I(lat_I<-0.416) = -0.416;
        lon_I = selected_site_lla(2)/180+sind(Az_selected_prn).*phi_E./cosd(lat_I*180); % semicircle

        % Get local time for tow_sequence
        local_time = mod(43200*lon_I+tow_sequence{epoch,selected_prn}, 3600*24);

        % Get IPP geomagnetic latitude
        lat_I_m = lat_I + 0.064*cosd((lon_I-1.617)*180);

        % Localize ionospheric parameters
        ionoparams = iono_all{epoch};

        if isempty(ionoparams)
            I_s{epoch,selected_prn} = NaN;
            continue;
        end

        % Klobuchar
        A1 = 5e-9;
        A3 = 50400;
        A2 = ionoparams(1)+ionoparams(2)*lat_I_m+ionoparams(3)*lat_I_m.^2+ionoparams(4)*lat_I_m.^3;
        A4 = ionoparams(5)+ionoparams(6)*lat_I_m+ionoparams(7)*lat_I_m.^2+ionoparams(8)*lat_I_m.^3;

        % Index for when SV is "in the dark" from the sun, appx
        index1 = abs(local_time-A3)<A4/4;

        % Sunlit epoch time steps
        I_z(index1) = (A1+A2(index1).*cos(2*pi*(local_time(index1)-A3)./A4(index1)))*c;
        % Sun-dark epoch time step
        % I_z(~index1) = A1*c;
        I_z(~index1) = NaN;

        % Rejecting non-visible PRN
        I_z(El_selected_prn<5) = nan;

        % Mapping function
        OF = 1 + 16*(0.53-El_selected_prn/180).^3;
        I_s{epoch,selected_prn} = I_z.*OF';

        % % Dual frequency (RINEX)
        % sTEC = f1^2*f2^2/40.3/(f1^2-f2^2)*(-C1C+C2L);
        % I_dual = 40.3*sTEC/f1^2;

    end
end

%%
%
% Now we have an 80x32 cell array of Klobuchar-modeled ionospheric delay
% vectors for each given PRN 1-32 and quarterly epoch from 2004 to 2023.
%
% Items where data was unavailable were left 0 for now.


% Visualize some entries
figure;
subplot(2,2,1)
hold on;
grid on;
for i = 1:80
    plot(linspace(0,24,96),I_s{i,8},'.-')
end
title(['PRN 8: Iono Delay, ' num2str(yyyy(1)) '-' num2str(yyyy(end))])
xlabel('GPS Time of Day')
ylabel('Ionospheric Correction [m]')


subplot(2,2,2)
hold on;
grid on;
for i = 1:28
    plot(linspace(0,24,96),I_s{i,8},'.-')
end
title(['PRN 8: Iono Delay, ' num2str(yyyy(1)) '-' num2str(yyyy(7))])
xlabel('GPS Time of Day')
ylabel('Ionospheric Correction [m]')

subplot(2,2,3)
hold on;
grid on;
for i = 29:52
    plot(linspace(0,24,96),I_s{i,8},'.-')
end
title(['PRN 8: Iono Delay, ' num2str(yyyy(8)) '-' num2str(yyyy(14))])
xlabel('GPS Time of Day')
ylabel('Ionospheric Correction [m]')

subplot(2,2,4)
hold on;
grid on;
for i = 53:80
    plot(linspace(0,24,96),I_s{i,8},'.-')
end
title(['PRN 8: Iono Delay, ' num2str(yyyy(15)) '-' num2str(yyyy(end))])
xlabel('GPS Time of Day')
ylabel('Ionospheric Correction [m]')

%%
% These plots show a high level overview of three distinct time periods of
% the whole, indicating visually that flux clearly varies with time.
%
% Converting these 1xN arrays into a surf plot, it is possible to get a 3D
% representation of these delays per epoch set.

% surf plot
I_s_surfable = [];
for j = 1:80
    if size(I_s{j,8}) <= 0
        I_s_surfable(:,j) = NaN;
    else
        I_s_surfable(:,j) = I_s{j,8};
    end
end


% 3D plotting of a specific PRN over time, 3D POV
figure
hold on
grid on
surf(year_sequence,linspace(0,24,96),I_s_surfable)
xlabel('Epoch [solstice/equinox]');
ylabel('Timeseries for Epoch [UTC Hours]');
zlabel('Ionospheric Delay [m]')
title('(Isometric) PRN 8: Ionospheric Delay, Min/Max over 2004-2023')
view([30 30 30])


% 3D plotting of a specific PRN over time, timescale POV
figure
hold on
grid on
surf(year_sequence,linspace(0,24,96),I_s_surfable)
xlabel('Epoch [solstice/equinox]');
ylabel('Timeseries for Epoch [UTC Hours]');
zlabel('Ionospheric Delay [m]')
title('(Timescale) PRN 8: Ionospheric Delay, Min/Max over 2004-2023')
view([0 -90 0])


%%
%
% From this angle, I can see the minima and maxima of the modeled
% ionospheric delay to GPS signals for 80 epochs, being 2 solstice and
% equinox pairs per year from 2004 to 2023.
%
% Looking at a few max and min values for each epoch, the median, and the
% mean, I can try to correlate patterns over time.  Standard deviation is
% also considered to demonstrate statistical significant of variance.

year_sequence = linspace(2004,2023,80);

I_s_mean = nan(80,1);
I_s_median = nan(80,1);
I_s_max = nan(80,1);
I_s_max5_mean = nan(80,1);
I_s_min = nan(80,1);
I_s_min5_mean = nan(80,1);
I_s_var = nan(80,1);
I_s_std = nan(80,1);
test_nan = nan(1,96);

for j = 1:80
    I = I_s{j,8};
    if length(I) > 1 && ~all(isnan(I))
        I_s_mean(j) = mean(I(~isnan(I)));
        I_s_median(j) = median(I(~isnan(I)));
        I_s_max(j) = max(I(~isnan(I)));
        I_s_max5_mean(j) = mean(maxk(I(~isnan(I)),5));
        I_s_min(j) = min(I(~isnan(I)));
        I_s_min5_mean(j) = mean(mink(I(~isnan(I)),5));
        I_s_var(j) = var(I(~isnan(I)));
        I_s_std(j) = std(I(~isnan(I)));
    else
        continue;
    end
end

figure
hold on;
grid on;
plot(year_sequence,I_s_mean)
plot(year_sequence,fillmissing(I_s_mean,'linear'))
plot(year_sequence,I_s_median)
plot(year_sequence,fillmissing(I_s_median,'linear'))
plot(year_sequence,I_s_max)
plot(year_sequence,fillmissing(I_s_max,'linear'))
plot(year_sequence,I_s_max5_mean)
plot(year_sequence,fillmissing(I_s_max5_mean,'linear'))
plot(year_sequence,I_s_min)
plot(year_sequence,fillmissing(I_s_min,'linear'))
plot(year_sequence,I_s_min5_mean)
plot(year_sequence,fillmissing(I_s_min5_mean,'linear'))
xlabel('Epoch [fractional year]')
ylabel('Ionospheric Delay [m]')
title('Mass Plot of Ionospheric Delay Characteristics')

%%
%
% Now, more of the story is being told.  It seems, before 2014, NASA CDDIS
% was publishing at a lower cadence.  Although we are smoothing over what
% is likely sub-cycles, calling |fillmissing| can help visualize trends
% with gaps in data.  It must be understood these missing terms are
% linearly interpolated, and capture the primary periodicity.


figure
hold on;
grid on;
plot(year_sequence,fillmissing(I_s_mean,'linear'),'.-')
plot(year_sequence,fillmissing(I_s_median,'linear'),'.-')
plot(year_sequence,fillmissing(I_s_max,'linear'),'.-')
plot(year_sequence,fillmissing(I_s_min,'linear'),'.-')
legend('Mean','Median','Max','Min')
xlabel('Epoch [fractional year]')
ylabel('Ionospheric Delay [m]')

%%
%
% Now, let's look at some variables in isolation directly.


% Fillmissing with linear interpolation
I_s_var_filled = fillmissing(I_s_var,"linear");
I_s_std_filled = fillmissing(I_s_std,"linear");
I_s_mean_filled = fillmissing(I_s_mean,'linear');

figure;
hold on;
grid on;
plot(year_sequence,I_s_var_filled)
xlabel('Epoch [fractional year]')
ylabel('Iono Delay Variance')
title('Variance of Ionospheric Delay over Time')


figure;
hold on;
grid on;
plot(year_sequence,I_s_mean_filled-I_s_std_filled,'--')
plot(year_sequence,I_s_mean_filled,'.-')
plot(year_sequence,I_s_mean_filled+I_s_std_filled,'--')
legend('\mu + \sigma','\mu','\mu - \sigma')
xlabel('Epoch [fractional year]')
ylabel('Iono Delay Variance')
title('Standard Deviation of Ionospheric Delay over Time')

%%
% This plot clearly confirms the variance observations.  The standard
% deviation is solved for 95% confidence of a presumed double-tailed normal
% distribution, so observations above that confidence interval indicate
% statistical anomalies of significance, and they can not be presumed
% intrinsic to the data itself and are likely induced by external
% phenomena.
%
% *This investigation serves to determine if that external influence is
% solar flux, and through these periodograms we are observing positive
% correlation.*


%%
%
% Periodicity of variance?

% Periodogram
figure;
periodogram(I_s_var_filled)

% Manual representation
abs_fft_iono_var = abs(fft(I_s_var_filled));
figure;
hold on;
grid on;
plot(abs_fft_iono_var(1:end/2))
title('Manual FFT of Ionospheric Delay Variance')
xlabel('Frequency')
ylabel('Ampitude')

%%
%
% These periodicity plots of the ionospheric delay variance over time are
% interesting.  They indicate there is a definite pattern which peaks
% around 0.5 times the epoch domain.  Since we are looking from 2004-2023,
% that would be a period of roughly *10 years*.
%
% In order to determine relevance of this periodicity against existing 
% solar data, a cross-correlation can be conducted between them to
% determine similarity.  According to my hypothesis, *years with greater
% variance in ionospheric delay will positively correlate with years
% experiencing higher than average solar flux.*

% [r, lag] = xcorr(I_s_var_filled,fluxTable.fluxobsflux);
[r, lag] = xcorr(fillmissing(I_s_max5_mean,'linear'),fluxTable.fluxobsflux);

figure
hold on;
grid on;
plot(linspace(0,19,length(r(1:ceil(end/2+1)))),r(1:ceil(end/2+1)))
xlabel('Shifted Years')
ylabel('Cross-Correlation Magnitude')
title('')

figure
hold on
grid on
periodogram(r(1:ceil(end/2+1)))

%%
%
% Both the cross-correlation and periodogram indicate that the solar flux
% intensity over time and ionospheric signal delay share a common frequency
% of half the time domain, or ~10 years.
%
% The cross-correlation shows how similar the trends are without any
% shifting near zero, and how similar they are again after 10 years of
% shifting, with an increasing trend after 1 years of shifting, indicating
% the next 2 years in the future from the existing data set may experience
% higher than average solar flux intensity.


%% Reference Signal Comparison

% Generate an 11-Year Period Signal

% Parameters
fs = 100;                % Sampling frequency (samples per year)
years = 0:1/fs:22;       % Time vector (22 years for 2 cycles)
period = 11;             % Period of the signal (years)
amplitude = 1;           % Amplitude of the signal
phase = 0;               % Phase offset (radians)
noise_level = 0.1;       % Noise level (amplitude)

% Generate the signal
signal = amplitude * sin(2 * pi * (1/period) * years + phase);

% Add noise to the signal
noisy_signal = signal + noise_level * randn(size(signal));

% Plot the signal
figure;
hold on; grid on;
plot(years, signal, 'b-', 'LineWidth', 2);  % Original signal
plot(years, noisy_signal, 'r.', 'MarkerSize', 5);  % Noisy signal
xlabel('Years');
ylabel('Signal Amplitude');
title('Synthetic 11-Year Period Signal');
legend('Original Signal', 'Noisy Signal');


%%
%
% Here, we are doing a direct cross correlation between observed flux and
% the control generated signal.  The cross correlation is then high-pass
% filtered to allow the rest of the spectrum to be observable

[r_ref, lag_ref] = xcorr(signal, fluxTable.fluxobsflux);
[r_noise_ref, lag_noise_ref] = xcorr(noisy_signal, fluxTable.fluxobsflux);

% high-pass filters
ref_hp = highpass(r_ref,.01);
ref_noise_hp = highpass(r_noise_ref,.01);

ref_fft = abs(fft(ref_hp));
ref_noise_fft = abs(fft(ref_noise_hp));

figure
periodogram(r_ref)

figure
hold on;
grid on;
plot(linspace(0,19,length(r_ref(1:ceil(end/2)))),r_ref(1:ceil(end/2)))
xlabel('Year')
ylabel('Amplitude')
title('Cross Correlation of Control and Klobuchar Delay')

%%
% These plots show very similar trends with cross-correlation peaks
% directly at the timespan cadence.  With the 11-year reference signal over
% a 19 year time span, the peak frequency from the periodogram and
% downselected FFT are half the total period, adding confidence to signal
% similarity.  Further, the variance increasing beyond standard deviation
% confirms the change in ionopsheric delay is statistically significant 
% when compared to solar flux.


%% Cleanup

fclose all;