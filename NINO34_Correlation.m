% Find ENSO-Sensitive records
% C.A. Morris, 2025

%% Load CoralHydro2k Database
% https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/noaa-coral-35453/html
load('CoralHydro2k_Matlab.mat');

%% Load SST data (ERSSTv5)
% ERSSTv5
% https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/
% Subset temporal coverage: 1854-01-01 to present
ERSST = 'ersst.mnmean.nc';
%ERSST = 'ERSSTv5_mnmean.nc';
ERSST_lat = ncread(ERSST,'lat');
ERSST_lon = ncread(ERSST,'lon');
ERSST_time = ncread(ERSST,'time');  % in days since reference
ERSST_sst = ncread(ERSST,'sst');
% Convert ERSST time to MATLAB datenum.
time_origin_ersst = datenum(1800,1,1);
ERSST_datenum = time_origin_ersst + ERSST_time;
ERSST_datetime = datetime(ERSST_datenum, 'ConvertFrom', 'datenum');
ERSST_datetime.Format = 'MM/dd/yyyy';

%% Filter Coral Data
% Remove year-only records
mask = strcmp({TS.paleoData_variableName}, 'year');
searchTS = TS(~mask);
% Remove d18O_sw records
mask = ~strcmp({searchTS.paleoData_variableName}, 'd18O_sw');
searchTS = searchTS(mask);
% Filter for Pacific Ocean records
mask = contains({searchTS.geo_ocean}, 'pacific', 'ignoreCase', true);
searchTS = searchTS(mask);
% Filter for Porites species
mask = contains({searchTS.paleoData_archiveSpecies}, 'Porites', 'ignoreCase', true);
searchTS = searchTS(mask);
% Filter for sub-annual resolution (exclude annual records)
mask_resolution = ~contains({searchTS.hasResolution_nominal}, 'annual', 'ignoreCase', true);
filtered_coral = searchTS(mask_resolution);
% Separate d18O and Sr/Ca records
mask_d18O = contains({filtered_coral.paleoData_variableName}, 'd18O', 'ignoreCase', true);
CH2K_d18O = filtered_coral(mask_d18O);

mask_srca = contains({filtered_coral.paleoData_variableName}, 'srca', 'ignoreCase', true);
CH2K_srca = filtered_coral(mask_srca);

% Filter for records extending past 1950
maskValid_d18O = arrayfun(@(x) max(x.year) >= 1950, CH2K_d18O);
CH2K_d18O = CH2K_d18O(maskValid_d18O);

maskValid_srca = arrayfun(@(x) max(x.year) >= 1950, CH2K_srca);
CH2K_srca = CH2K_srca(maskValid_srca);

%% Compute coral anomalies (remove seasonal cycle) 
function [dt, anom] = Coral_Anom(coral_time, coral_vals)
    % Convert fractional years to year and month
    years = floor(coral_time);
    months = floor((coral_time - years) * 12) + 1;
    months(months < 1) = 1;
    months(months > 12) = 12;
    % Compute monthly climatology
    seasonal_mean = NaN(12,1);
    for m = 1:12
        idx = (months == m);
        if any(idx)
            seasonal_mean(m) = mean(coral_vals(idx), 'omitnan');
        end
    end
    % Remove seasonal cycle
    coral_anom = coral_vals - seasonal_mean(months);
    % Create datetime vector (assign day 15 of each month)
    coral_dt = datetime(years, months, 15);
    [dt, sortIdx] = sort(coral_dt);
    anom = coral_anom(sortIdx);
end

%% Process coral records
% d18O records
for i = 1:length(CH2K_d18O)
    [CH2K_d18O(i).dt, CH2K_d18O(i).anom] = Coral_Anom(...
        CH2K_d18O(i).year, CH2K_d18O(i).paleoData_values);
end

% Sr/Ca records
for i = 1:length(CH2K_srca)
    [CH2K_srca(i).dt, CH2K_srca(i).anom] = Coral_Anom(...
        CH2K_srca(i).year, CH2K_srca(i).paleoData_values);
end

%% Compute NINO3.4 SST anomalies
time_mask = ERSST_datetime >= datetime(1950,1,1);
ERSST_datetime = ERSST_datetime(time_mask);
ERSST_sst = ERSST_sst(:,:,time_mask);

% Define NINO3.4 region (longitude in degrees East)
% NINO3.4: 5°N-5°S, 170°W-120°W (190°E-240°E)
nino34_lat = [-5, 5];
nino34_lon = [190, 240];

% Extract NINO3.4 ERSST
lat_idx = find(ERSST_lat >= nino34_lat(1) & ERSST_lat <= nino34_lat(2));
lon_idx = find(ERSST_lon >= nino34_lon(1) & ERSST_lon <= nino34_lon(2));
nino34_sst = ERSST_sst(lon_idx, lat_idx, :);
% Compute spatial average over NINO3.4 region
nino34_sst_ts = squeeze(mean(mean(nino34_sst, 1, 'omitnan'), 2, 'omitnan'));

% Compute monthly climatology
months = month(ERSST_datetime);
nino34_climatology = NaN(12, 1);
for m = 1:12
    month_idx = (months == m);
    nino34_climatology(m) = mean(nino34_sst_ts(month_idx), 'omitnan');
end

% Remove seasonal cycle
nino34_ssta = NaN(size(nino34_sst_ts));
for i = 1:length(nino34_sst_ts)
    nino34_ssta(i) = nino34_sst_ts(i) - nino34_climatology(months(i));
end

%% Correlation function w/ effective DOF
% Account for autocorrelation
function [r, p, dof_eff] = corr_eff(data1, data2, L)
    % Calculate correlation w/ effective DOF
    % Inputs:
    % data1: coral timeseries
    % data2: SST timeseries
    % L: lag cutoff for autocorrelation (default=1)
    % Outputs: 
    % r: correlation coefficient
    % p: significance level
    % dof_eff: effective degrees of freedom

    if nargin < 3
    L = 1;
    end
   
    % Put data in column vector
    data1 = data1(:);
    data2 = data2(:);
    valid = ~isnan(data1) & ~isnan(data2);
    data1 = data1(valid);
    data2 = data2(valid);
    N = length(data1);
    
    % Check minimum sample size
    if N < 3
        warning('Too few data points');
        r = NaN; p = NaN; dof_eff = NaN;
        return;
    end
    
    % Calculate correlation
    r = corr(data1, data2);
    % Get autocorrelations up to lag L
    r1 = zeros(1, L);
    r2 = zeros(1, L);
    
    for k = 1:L
        if N > k
            % Lag-k autocorrelation for data1
            c1 = corrcoef(data1(1:N-k), data1(k+1:N));
            r1(k) = c1(1,2);
            % Lag-k autocorrelation for data2
            c2 = corrcoef(data2(1:N-k), data2(k+1:N));
            r2(k) = c2(1,2);
        end
    end
    
    % Average autocorrelation
    r1_mean = mean(r1);
    r2_mean = mean(r2);
    
    % Calculate effective DOF for each series (Bretherton et al., 1999)
    dof_eff1 = N * (1 - r1_mean) / (1 + r1_mean);
    dof_eff2 = N * (1 - r2_mean) / (1 + r2_mean);
    
    % Use geometric mean of the two effective DOFs
    dof_eff = sqrt(dof_eff1 * dof_eff2);
    
    % Ensure minimum DOF for t-test
    dof_eff = max(dof_eff, 3);
    
    % Calculate p-values with t-test
    % Degrees of freedom for t-distribution
    df = dof_eff - 2;
    
    % Calculate t-statistic
    t_stat = abs(r) * sqrt(df) / sqrt(1 - r^2);
    
    % Calculate two-tailed p-value
    p = 2 * tcdf(-abs(t_stat), df);

end

%% Find correlations b/w coral anomalies & NINO3.4 ssta
results = table();
L = 1; % Adjust based on autocorrelation of data

% process d18O records
for i = 1:length(CH2K_d18O)
    % Get coral time series
    dt_coral = CH2K_d18O(i).dt;
    anom_coral = CH2K_d18O(i).anom;
    % Filter for 1950 onwards
    time_mask = dt_coral >= datetime(1950,1,1);
    dt_coral = dt_coral(time_mask);
    anom_coral = anom_coral(time_mask);
    % Interpolate SST anomalies to coral resolution
    sst34_interp = interp1(ERSST_datetime, nino34_ssta, dt_coral, 'linear', NaN);
    % Compute correlations
    [r34, p34, dof34] = corr_eff(anom_coral, sst34_interp, L);
    % Get temporal coverage
    coverageStart = min(CH2K_d18O(i).year);
    coverageEnd = max(CH2K_d18O(i).year);
    n_points = sum(~isnan(anom_coral) & ~isnan(sst34_interp));
    
    % Create results row for each record
    row = table(...
        {CH2K_d18O(i).dataSetName}', ...
        {CH2K_d18O(i).geo_siteName}', ...
        CH2K_d18O(i).geo_latitude, ...
        CH2K_d18O(i).geo_longitude, ...
        {CH2K_d18O(i).hasResolution_nominal}', ...
        {CH2K_d18O(i).paleoData_variableName}', ...
        coverageStart, ...
        coverageEnd, ...
        n_points, ...
        r34, ...
        p34, ...
        dof34, ...
        'VariableNames', {'dataSetName', 'geo_siteName', 'geo_latitude', 'geo_longitude', ...
        'resolution', 'variableName', 'coverageStart', 'coverageEnd', ...
        'n_datapoints', 'r_nino34', 'p_nino34', 'dof_nino34'});

    results = [results; row];
end

% process Sr/Ca records
for i = 1:length(CH2K_srca)
    % Get coral time series
    dt_coral = CH2K_srca(i).dt;
    anom_coral = CH2K_srca(i).anom;
    % Filter for 1950 onwards
    time_mask = dt_coral >= datetime(1950,1,1);
    dt_coral = dt_coral(time_mask);
    anom_coral = anom_coral(time_mask);
    % Interpolate SST anomalies to coral resolution
    sst34_interp = interp1(ERSST_datetime, nino34_ssta, dt_coral, 'linear', NaN);
    % Compute correlations
    [r34, p34, dof34] = corr_eff(anom_coral, sst34_interp, L);
    % Get temporal coverage
    coverageStart = min(CH2K_srca(i).year);
    coverageEnd = max(CH2K_srca(i).year);
    n_points = sum(~isnan(anom_coral) & ~isnan(sst34_interp));
    
    % Create results row for each record
    row = table(...
        {CH2K_srca(i).dataSetName}', ...
        {CH2K_srca(i).geo_siteName}', ...
        CH2K_srca(i).geo_latitude, ...
        CH2K_srca(i).geo_longitude, ...
        {CH2K_srca(i).hasResolution_nominal}', ...
        {CH2K_srca(i).paleoData_variableName}', ...
        coverageStart, ...
        coverageEnd, ...
        n_points, ...
        r34, ...
        p34, ...
        dof34, ...
        'VariableNames', {'dataSetName', 'geo_siteName', 'geo_latitude', 'geo_longitude', ...
        'resolution', 'variableName', 'coverageStart', 'coverageEnd', ...
        'n_datapoints', 'r_nino34', 'p_nino34', 'dof_nino34'});
    
    results = [results; row];
end