% Calculate interannual variance from coral d18O, Sr/Ca and SST
% Using 20-yeard standard deviation of 2-7 year bandpassed data
% C.A. Morris, 2025

%% FFT Bandpass function
function y=bandpass(x,l,h,fs)
%function y=bandpass(x,l,h,fs) bandpasses time series x.
%all frequencies higher than lpf and lower than hpf in cycles per time step
%are removed.
% x = input time series
% l = lower period bound (years)
% h = upper period bound (years) 
% fs = sampling frequency (samples per year)
% first lower #, then higher
if nargin < 4
    fs = 12; % default to monthly
end
lpf=1/l/fs; % Low pass frequency (cycles per time step)
hpf=1/h/fs; % High pass frequency (cycles per time step)
[len,num]=size(x);
span=len-1;
slope=(x(len,:)-x(1,:))/span;
off=(x(1,:)*len-x(len,:))/span;
match=(1:len)'*slope+ones(len,1)*off;
xi=fft(x-match);
mlpf=ceil(lpf*len);
mhpf=floor(hpf*len);
xi(1:mhpf+1,:)=zeros(mhpf+1,num);
xi(len-mhpf+1:len,:)=zeros(mhpf,num);
xi(mlpf+1:len-mlpf+1,:)=zeros(size(xi(mlpf+1:len-mlpf+1,:)));
y=ifft(xi);
end

%% Put ENSO-Sensitive CH2K records into a structure
coralTables = struct( ...
  'Kiritimati_d18O', Kiritimati_d18O, ...
  'Palmyra_d18O', Palmyra_d18O, ...
  'Palmyra_srca', Palmyra_srca, ...
  'Bunaken_d18O', Bunaken_d18O, ...
  'Padang_d18O', Padang_d18O, ...
  'Sabine_d18O', Sabine_d18O, ...
  'Sabine_srca', Sabine_srca, ...
  'Rarotonga1_d18O', Rarotonga1_d18O, ...
  'Rarotonga1_srca', Rarotonga1_srca, ...
  'Vanua_d18O', Vanua_d18O, ...
  'Vanua_srca', Vanua_srca, ...
  'Rarotonga2_d18O', Rarotonga2_d18O,  ...
  'Ulong_d18O', Ulong_d18O, ...
  'Clarion_d18O', Clarion_d18O, ...
  'Madang_d18O', Madang_d18O, ...
  'Laing_d18O', Laing_d18O, ...
  'Maiana_d18O', Maiana_d18O, ...
  'Clipperton_d18O', Clipperton_d18O, ...
  'Clipperton_srca', Clipperton_srca, ...
  'Amedee_srca', Amedee_srca, ...
  'Spliced19_d18O', Spliced19_d18O, ...
  'Spliced20_d18O', Spliced20_d18O, ...
  'X12_FS13_11', X12_FS13_11, ...
  'X14_FS28_10', X14_FS28_10 ...
);

%% Load and filter SST data
ERSST = 'ersst.mnmean.nc';
lat = ncread(ERSST,'lat');
lon = ncread(ERSST,'lon');
t_sst = ncread(ERSST,'time');
sst_mtx = ncread(ERSST,'sst');
dt_sst = datetime(1800,1,1) + days(t_sst);
lat_idx = lat>=-5 & lat<=5;
lon_idx = lon>=190 & lon<=240;
mtx = sst_mtx(lon_idx,lat_idx,:);
mnth = month(dt_sst);

% Remove seasonal cycle
clim = zeros(size(mtx,1), size(mtx,2), 12);
for m = 1:12
    clim(:,:,m) = mean(mtx(:,:,mnth==m), 3, 'omitnan');
end
N = size(mtx,3);
ERSST_ano = nan(size(mtx));
for i = 1:N
    ERSST_ano(:,:,i) = mtx(:,:,i) - clim(:,:,mnth(i));
end

% Niño3.4 anomalies
nino34_raw = squeeze(mean(mean(ERSST_ano,1,'omitnan'),2,'omitnan'));

% 10-year high-pass filter
nino34_decadal = movmean(nino34_raw, 120, 'omitnan');
nino34_hp = nino34_raw - nino34_decadal;

% 2–7 yr FFT bandpass
nino34_bp = bandpass(nino34_hp, 2, 7, 12); % fs=12 for monthly SST data
nino34_bp = real(nino34_bp);

% 20 yr running standard deviation
nSD = movstd(nino34_bp, 20*12, 'omitnan', 'Endpoints', 'fill');

% Caculate Baseline
idxS = year(dt_sst) >= 1969 & year(dt_sst) <= 1989;
baseS = mean(nSD(idxS), 'omitnan');

% Calculate % change in SST variance
nino_pct = ((nSD - baseS) / baseS) * 100;

%% Process Coral Data
% Parameters
baseline_start = 1969;
baseline_end = 1989;
filter_order = 2;
low_cut = 1/7;
high_cut = 1/2;
win_yrs = 20;

records = { ...
  'Kiritimati_d18O','Palmyra_d18O','Palmyra_srca','Bunaken_d18O',...
  'Padang_d18O','Sabine_d18O','Sabine_srca','Rarotonga1_d18O',...
  'Rarotonga1_srca','Vanua_d18O','Vanua_srca','Rarotonga2_d18O',...
  'Ulong_d18O','Clarion_d18O','Madang_d18O','Laing_d18O',...
  'Maiana_d18O','Clipperton_d18O','Clipperton_srca',...
  'Amedee_srca','Spliced19_d18O','Spliced20_d18O',...
  'X12_FS13_11','X14_FS28_10' ...
};

% Account for different data resolution
pts_per_year = containers.Map();
monthly = { ...
  'Kiritimati_d18O','Palmyra_d18O','Palmyra_srca','Bunaken_d18O',...
  'Padang_d18O','Sabine_d18O','Sabine_srca','Ulong_d18O',...
  'Amedee_srca','Clipperton_d18O','Clipperton_srca','Spliced19_d18O','Spliced20_d18O' ...
};
for m = monthly
    pts_per_year(m{1}) = 12;
end
pts_per_year('Rarotonga1_d18O') = 9;
pts_per_year('Rarotonga1_srca') = 9;
pts_per_year('Vanua_d18O') = 8;
pts_per_year('Vanua_srca') = 8;
pts_per_year('Rarotonga2_d18O') = 8;
pts_per_year('Maiana_d18O') = 6;
for nm = ["Clarion_d18O","Madang_d18O","Laing_d18O"]
    pts_per_year(char(nm)) = 4;
end
% Uneven fossil corals
yrs12 = max(X12_FS13_11{:,1}) - min(X12_FS13_11{:,1});
pts_per_year('X12_FS13_11') = round(size(X12_FS13_11,1)/yrs12);
yrs14 = max(X14_FS28_10{:,1}) - min(X14_FS28_10{:,1});
pts_per_year('X14_FS28_10') = round(size(X14_FS28_10,1)/yrs14);

baselineRecs = {'Spliced19_d18O','Spliced20_d18O','X12_FS13_11','X14_FS28_10'};
singleRecs   = {'Spliced20_d18O','X12_FS13_11','X14_FS28_10'};

coral = struct();
for k = 1:numel(records)
    rec = records{k};
    T   = coralTables.(rec);
    tf  = T{:,1}; % fractional years
    obs = T{:,2}; % raw d18O or Sr/Ca
    fs  = pts_per_year(rec);
    % convert to datetime
    dt  = datetime(floor(tf),1,1) + days(round((tf-floor(tf))*365.25));
    
    coral.(rec).tf = tf;
    coral.(rec).time = dt;
    coral.(rec).fs = fs;
    coral.(rec).d18O_13mo = movmean(obs, 13, 'omitnan');
end

%% Calculate Baseline
x_tf  = coral.Kiritimati_d18O.tf;
x_obs = coralTables.Kiritimati_d18O{:,2};
x_fs  = coral.Kiritimati_d18O.fs;
x_dt  = coral.Kiritimati_d18O.time;

x_month = month(x_dt);
x_ano = nan(size(x_obs));
for m = 1:12
    idx = x_month == m;
    x_ano(idx) = x_obs(idx) - mean(x_obs(idx), 'omitnan');
end
x_decadal = movmean(x_ano, 120, 'omitnan');
x_hp = x_ano - x_decadal;

% FFT Bandpass
x_bp = bandpass(x_hp, 2, 7, x_fs);
x_bp = real(x_bp);
x_std20 = movstd(x_bp, win_yrs * x_fs, 'omitnan');

idx_baseline = x_tf >= baseline_start & x_tf <= baseline_end;
globalBaseline_std = mean(x_std20(idx_baseline), 'omitnan');

%% Calculate Uncertainty for short xmas fossil coral
variance_err_sym = zeros(numel(singleRecs),1);
x13mo = coral.Kiritimati_d18O.d18O_13mo;
x_tf  = coral.Kiritimati_d18O.tf;
idx_baseline_short = x_tf >= baseline_start & x_tf <= baseline_end;
baseline_std_short = std(x13mo(idx_baseline_short), 'omitnan');

for i = 1:numel(singleRecs)
    rec = singleRecs{i};
    tf_rec = coral.(rec).tf;
    tf_rec = tf_rec(~isnan(tf_rec));
    rec_yrs = max(tf_rec) - min(tf_rec);
    T_yrs = round(rec_yrs);
    Tpts = max(T_yrs*12,1);

    std_T = movstd(x13mo, Tpts, 'omitnan');
    std_20 = movstd(x13mo, 240,   'omitnan');
    valid_idx = ~isnan(std_T) & ~isnan(std_20);
    diff_std = abs(std_T(valid_idx) - std_20(valid_idx));
    var_uncert = mean(diff_std);

    variance_err_sym(i) = (var_uncert / baseline_std_short) * 100;
end

%% Filter all records
for k = 1:numel(records)
    rec = records{k};
    tf  = coral.(rec).tf;
    obs = coralTables.(rec){:,2};
    fs  = coral.(rec).fs;
    dt  = coral.(rec).time;

    valid = ~isnan(tf) & ~isnan(obs);
    tf  = tf(valid);
    obs = obs(valid);
    dt  = dt(valid);

    record_len_yrs = max(tf) - min(tf);

    if ismember(rec, singleRecs) || record_len_yrs < 20
        % For records <20 years take 13 month running avg
        sig = std(coral.(rec).d18O_13mo, 'omitnan');
        pct = (sig - baseline_std_short) / baseline_std_short * 100;

        coral.(rec).single_pct = pct;
        coral.(rec).mid_time = mean(tf,'omitnan');
        idx_s = find(strcmp(singleRecs, rec));
        v = variance_err_sym(idx_s);
        coral.(rec).single_err = [v, v];
        coral.(rec).std_run = sig * ones(size(tf));
        coral.(rec).pct_change = pct * ones(size(tf));
        continue
    end

    mVec = month(dt);
    d18O_ano = nan(size(obs));
    for mm = 1:12
        idx_mm = (mVec == mm);
        d18O_ano(idx_mm) = obs(idx_mm) - mean(obs(idx_mm), 'omitnan');
    end
    % 10 year highpass
    d18O_decadal = movmean(d18O_ano, round(10*fs), 'omitnan');
    d18O_hp      = d18O_ano - d18O_decadal;

    % 2-7 year bandpass
    filt = bandpass(d18O_hp, 2, 7, fs);
    filt = real(filt);
    run_sd = movstd(filt, win_yrs*fs, 'omitnan');

    coral.(rec).filt    = filt;
    coral.(rec).std_run = run_sd;

    % Basleline standard deviation for each site
    if ismember(rec, baselineRecs)
        coral.(rec).base_std = globalBaseline_std;
    else
        idx_b = tf >= baseline_start & tf <= baseline_end;
        coral.(rec).base_std = mean(run_sd(idx_b), 'omitnan');
    end

    % Calcualte % change in standard deviation
    coral.(rec).pct_change = (run_sd - coral.(rec).base_std) ./ coral.(rec).base_std * 100;
end

%% Interpolate all records to a common bimonthly time grid
tstart = datetime(1815,1,1);
tend = datetime(2025,12,1);
Tm = (tstart:calmonths(2):tend)';

plotFields = {'Kiritimati_d18O','Spliced19_d18O','Maiana_d18O','Padang_d18O','Amedee_srca'};
useRecs = [plotFields singleRecs];
R = numel(useRecs);

Mstd = nan(numel(Tm), R);
Mpct = nan(numel(Tm), R);

for k = 1:R
    rec = useRecs{k};
    if ismember(rec, singleRecs)
        % single‐point fossil
        mid_tf = coral.(rec).mid_time;
        fracTm = Tm.Year + (Tm.Month-0.5)/12;
        [~, idx0] = min(abs(fracTm - mid_tf));
        Mstd(idx0,k) = coral.(rec).std_run(1);
        Mpct(idx0,k) = coral.(rec).single_pct;
    else
        % Interpolate onto bimonthly grid
        tfFull = coral.(rec).tf;
        Ystd = coral.(rec).std_run;
        Ypct = coral.(rec).pct_change;
        [Xu, ia] = unique(tfFull);
        Mstd(:,k) = interp1(Xu, Ystd(ia), Tm.Year+(Tm.Month-0.5)/12, 'linear', NaN);
        Mpct(:,k) = interp1(Xu, Ypct(ia), Tm.Year+(Tm.Month-0.5)/12, 'linear', NaN);
    end
end

% NINO3.4 SST
Xsst = dt_sst.Year + (dt_sst.Month-0.5)/12;
[Xu_s, ia_s] = unique(Xsst);
nino_pct6 = interp1(Xu_s, nino_pct(ia_s), Tm.Year+(Tm.Month-0.5)/12, 'linear', NaN);

fullPctCols = 1:numel(plotFields);
MpctFull = Mpct(:, fullPctCols);
crossMeanPct_tm = mean(MpctFull, 2, 'omitnan');
numPctRec_tm = sum(~isnan(MpctFull), 2);
crossSEPct_tm  = std(MpctFull, 0, 2, 'omitnan') ./ sqrt(numPctRec_tm);
crossMeanPct_run = movmean(crossMeanPct_tm, 6, 'omitnan');
crossSEPct_run = movmean(crossSEPct_tm,   6, 'omitnan');

%% Plot
figure('Color','w','Units','inches','Position',[1 1 9 7]);
set(gcf, ...
    'DefaultAxesFontName','Helvetica', ...
    'DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontSize',16, ...
    'DefaultTextFontSize',16, ...
    'DefaultLegendInterpreter','tex' ...
);
hold on;

% 0 baseline
yline(0, '--', ...
    'Color', [0.3 0.3 0.3], ...
    'LineWidth', 1.5, ...
    'HandleVisibility','off');

% +/- 2 SE envelope (excluding short fossil corals)
envMask = ~isnan(crossMeanPct_run) & ~isnan(crossSEPct_run);
Tm_env = Tm(envMask);
upperEnv = crossMeanPct_run(envMask) + 2*crossSEPct_run(envMask);
lowerEnv = crossMeanPct_run(envMask) - 2*crossSEPct_run(envMask);

patch( ...
    [Tm_env; flipud(Tm_env)], ...
    [upperEnv; flipud(lowerEnv)], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.5, ...
    'HandleVisibility','off' ...
);

%
legendNames = { ...
    'Kiritimati \delta^{18}O', ... % xmas d18O
    '',                             ... % Spliced19_d18O 
    'Maiana, Kiribati \delta^{18}O', ...
    'Padangbai, Bali \delta^{18}O', ...
    'Amédée Island Sr/Ca'        ...
};

% Colors
cols = { ...
    [0 0 1], ...
    [0 0 1], ...
    [0.9, 0.2, 0.6], ...
    [0.85, 0.33, 0.1], ...
    [0.6, 0.3, 0.7] ...
};

% Plot all coral % change
for i = 3:numel(plotFields)
    if isempty(legendNames{i})
        plot(Tm, Mpct(:,i), '-', ...
            'Color', cols{i}, ...
            'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    else
        plot(Tm, Mpct(:,i), '-', ...
            'Color', cols{i}, ...
            'LineWidth', 1.5, ...
            'DisplayName', legendNames{i});
    end
end

for i = 1:2
    if isempty(legendNames{i})
        plot(Tm, Mpct(:,i), '-', ...
            'Color', cols{i}, ...
            'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    else
        plot(Tm, Mpct(:,i), '-', ...
            'Color', cols{i}, ...
            'LineWidth', 2.5, ...
            'DisplayName', legendNames{i});
    end
end

% Fossil coral errorbars
for ii = 1:numel(singleRecs)
    recPt = singleRecs{ii};
    fracTm = Tm.Year + (Tm.Month-0.5)/12;
    [~, idx0] = min(abs(fracTm - coral.(recPt).mid_time));
    yVal = coral.(recPt).single_pct;
    e    = coral.(recPt).single_err;
    errorbar(Tm(idx0), yVal, e(1), e(2), 'o', ...
             'LineStyle', 'none', 'LineWidth', 2, 'MarkerSize', 10, ...
             'Color', [0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], ...
             'HandleVisibility', 'off');
end

% Nino 3.4 SST % change
plot(Tm, nino_pct6, '-', ...
     'Color', [0.4 0.4 0.4], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'NIÑO3.4 SST (ERSSTv5)');

% Coral running average
numOverlap = sum(~isnan(MpctFull), 2);
idx_run    = numOverlap > 1;
plot(Tm(idx_run), crossMeanPct_run(idx_run), 'k-', ...
     'LineWidth', 2, ...
     'DisplayName', 'Running Average Variance');


xticks_vec = datetime(1820:20:2025, 1, 1);

% Axes formatting
xlim([tstart, datetime(2025,12,1)]);
ylim([-80, 40]);
set(gca, ...
    'FontName', 'Helvetica', ...
    'FontSize', 16, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'XTick', xticks_vec, ...
    'YTick', -80:20:80 ...
);
xtickangle(0);

% Axis labels
xlabel('Year (CE)', 'FontName', 'Helvetica', 'FontSize', 18);
ylabel('Percent change in variance (%)', 'FontName', 'Helvetica', 'FontSize', 18);

% Legend
lg = legend('Orientation', 'horizontal', ...
            'NumColumns', 2, ...
            'Location', 'northoutside');
set(lg, 'FontName', 'Helvetica', 'FontSize', 16, 'Box', 'on');

box on