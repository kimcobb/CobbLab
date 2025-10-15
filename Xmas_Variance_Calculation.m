% Calculate interannual variance from coral d18O and SST
% Using standard deviation of 13 month running average
% C.A. Morris, 2025

%% Parameters
baseline = [1987, 2007];
fs = 12; % monthly resolution
win_pts = 20 * fs; % 20-year window in points

%% Load coral d18O data files
coral_files = {
    'Christmas_Composited18O.csv', ...
    'SplicedRecord_19cent.csv', ...
    'SplicedRecord_20cent.csv', ...
    'X12_FS13_11.csv', ...
    'X14_FS28_10.csv'
};
composite_file = 'Christmas_Composited18O.csv';

%% Process coral records
records = struct();
baseline_std = NaN;

for i = 1:length(coral_files)
    filename = coral_files{i};
    data = readtable(filename);
    time = data.time;
    d18O = data.d18O;
    record_length = max(time) - min(time);
    
    % Convert fractional years to months
    year_part = floor(time);
    month_num = round((time - year_part) * 12 + 0.5);
    month_num(month_num < 1) = 1;
    month_num(month_num > 12) = 12;
    
    % Remove seasonal cycle
    d18O_ano = d18O;
    for m = 1:12
        idx = (month_num == m);
        if any(idx)
            d18O_ano(idx) = d18O(idx) - mean(d18O(idx), 'omitnan');
        end
    end
    
    % 10-year high-pass filter
    d18O_decadal = movmean(d18O_ano, 120, 'omitnan');
    d18O_hp = d18O_ano - d18O_decadal;
    
    % 13-month moving average
    d18O_filt = movmean(d18O_hp, 13, 'omitnan');
    
    % Calculate 20-year sliding standard deviation
    d18O_std20 = movstd(d18O_filt, win_pts, 'omitnan');
    
    % Get baseline standard deviation from composite record
    if strcmp(filename, composite_file)
        idx_baseline = time >= baseline(1) & time <= baseline(2);
        baseline_std = mean(d18O_std20(idx_baseline), 'omitnan');
    end
    
    % Calculate percent change from baseline
    pct_change = ((d18O_std20 - baseline_std) / baseline_std) * 100;
    
    % Store results
    record_name = strrep(filename, '.csv', '');
    record_name = strrep(record_name, '.', '_'); % Replace dots with underscores
    records.(record_name) = struct(...
        'time', time, ...
        'd18O', d18O, ...
        'd18O_ano', d18O_ano, ...
        'd18O_hp', d18O_hp, ...
        'd18O_filt', d18O_filt, ...
        'std20', d18O_std20, ...
        'pct_change', pct_change, ...
        'record_length', record_length ...
    );
end

%% Calculate uncertainty for short records using composite
% Use composite record to estimate uncertainty for short records
composite_name = strrep(composite_file, '.csv', '');
composite_name = strrep(composite_name, '.', '_');
comp = records.(composite_name);
comp_time = comp.time;
comp_d18O = comp.d18O;

% Convert fractional years to months
year_part = floor(comp_time);
month_num = round((comp_time - year_part) * 12 + 0.5);
month_num(month_num < 1) = 1;
month_num(month_num > 12) = 12;

% Remove seasonal cycle
comp_deseasoned = comp_d18O;
for m = 1:12
    idx = (month_num == m);
    if any(idx)
        comp_deseasoned(idx) = comp_d18O(idx) - mean(comp_d18O(idx), 'omitnan');
    end
end

% 10-year high-pass filter
comp_hp = comp_deseasoned - movmean(comp_deseasoned, 120, 'omitnan');

% 13-month moving average
comp_13 = movmean(comp_hp, 13, 'omitnan');

% Calculate uncertainty for each short record
variance_err_sym = struct();
record_names = fieldnames(records);

for i = 1:length(record_names)
    name = record_names{i};
    rec_length = records.(name).record_length;
    
    if rec_length < 20
        % Compare variance estimates using different window lengths
        T_yrs = round(rec_length);
        T_pts = max(T_yrs * fs, 1);
        std_T = movstd(comp_13, T_pts, 'omitnan');
        std_20 = movstd(comp_13, win_pts, 'omitnan');
        
        % Uncertainty = mean difference between estimates
        valid = ~isnan(std_T) & ~isnan(std_20);
        if any(valid)
            uncertainty = mean(abs(std_T(valid) - std_20(valid)));
            variance_err_sym.(name) = (uncertainty / baseline_std) * 100;
        else
            variance_err_sym.(name) = NaN;
        end
    end
end

%% Process SST data (ERSST)
sst_file = 'ersst.mnmean.nc';
sst_lat = ncread(sst_file, 'lat');
sst_lon = ncread(sst_file, 'lon');
sst_time_raw = ncread(sst_file, 'time');
sst = ncread(sst_file, 'sst');
    
% Convert time
dt_sst = datetime(1800, 1, 1) + days(sst_time_raw);
sst_time = year(dt_sst) + (month(dt_sst) - 0.5) / 12;
    
% Extract NINO3.4 region (5S-5N, 190E-240E)
lat_idx = sst_lat >= -5 & sst_lat <= 5;
lon_idx = sst_lon >= 190 & sst_lon <= 240;
sst_data = squeeze(mean(mean(sst(lon_idx, lat_idx, :), 1, 'omitnan'), 2, 'omitnan'));
    
% Remove seasonal cycle
month_sst = month(dt_sst);
sst_ano = sst_data;
for m = 1:12
    idx = (month_sst == m);
    if any(idx)
        sst_ano(idx) = sst_data(idx) - mean(sst_data(idx), 'omitnan');
    end
end
    
% 10 year High-pass filte
sst_hp = sst_ano - movmean(sst_ano, 120, 'omitnan');
    
% 13-month moving average
sst_bp = movmean(sst_hp, 13, 'omitnan');
    
% 20-year running standard deviation
sst_std20 = movstd(sst_bp, win_pts, 'omitnan');
    
% Baseline standard deviation
idx_baseline_sst = sst_time >= baseline(1) & sst_time <= baseline(2);
sst_baseline_std = mean(sst_std20(idx_baseline_sst), 'omitnan');
    
 % Percent change from baseline
sst_pct_change = ((sst_std20 - sst_baseline_std) / sst_baseline_std) * 100;


%% Figure setup
time_segments = [1820 1851; 1888 1900; 1912 2025];
n_segments = size(time_segments, 1);

record_colors = [
    0, 0, 0; % Christmas_Composited18O
    0, 0, 1; % SplicedRecord_19cent
    0, 0, 1; % SplicedRecord_20cent
    0, 0, 0; % X12-FS13-11
    0, 0, 0 % X14-FS28-10
];
sst_color = [0.5, 0.5, 0.5];

figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 11 5]);

% Layout settings
left_margin = 0.15;
right_margin = 0.15;
gap = 0.02;
seg_lengths = time_segments(:,2) - time_segments(:,1);
total_years = sum(seg_lengths);
total_width = 1 - left_margin - right_margin;
seg_widths = (total_width - gap*(n_segments-1)) * seg_lengths / total_years;

%% Plot by time segment
for seg = 1:n_segments
    if seg == 1
        left_pos = left_margin;
    else
        left_pos = left_margin + sum(seg_widths(1:seg-1)) + gap*(seg-1);
    end

    ax(seg) = axes('Position', [left_pos, 0.12, seg_widths(seg), 0.78]); 
    hold on;

    % 0 Reference line
    plot([time_segments(seg,1), time_segments(seg,2)], [0, 0], '--', ...
        'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.5);

    % Coral records
    for i = 1:length(record_names)
        name = record_names{i};
        t = records.(name).time;
        y = records.(name).pct_change;
        rec_length = records.(name).record_length;

        in_segment = t >= time_segments(seg,1) & t <= time_segments(seg,2);

        if rec_length >= 20
            % Long records as a continuous line
            valid_idx = in_segment & ~isnan(y);
            if any(valid_idx)
                plot(t(valid_idx), y(valid_idx), '-', ...
                    'Color', record_colors(i,:), 'LineWidth', 2);
            end
        else
            % Short records as a point with uncertainty bars
            if any(in_segment)
                mid_time = mean(t(in_segment));
                mean_pct = mean(y(in_segment), 'omitnan');

                % Pull uncertainty
                if exist('variance_err_sym', 'var') && ...
                   isfield(variance_err_sym, name) && ...
                   ~isnan(variance_err_sym.(name))

                    err = variance_err_sym.(name);
                    errorbar(mid_time, mean_pct, err, err, 'o', ...
                        'Color', record_colors(i,:), ...
                        'MarkerFaceColor', record_colors(i,:), ...
                        'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 8);
                else
                    plot(mid_time, mean_pct, 'o', ...
                        'Color', record_colors(i,:), ...
                        'MarkerFaceColor', record_colors(i,:), ...
                        'MarkerSize', 8);
                end
            end
        end
    end

    % SST
    sst_in_seg = sst_time >= time_segments(seg,1) & sst_time <= time_segments(seg,2);
    valid_sst = sst_in_seg & ~isnan(sst_pct_change);
    if any(valid_sst)
        plot(sst_time(valid_sst), sst_pct_change(valid_sst), '-', ...
            'Color', sst_color, 'LineWidth', 2);
    end

    % Axes formatting
    xlim(time_segments(seg,:));
    ylim([-80, 80]);
    set(gca, 'FontSize', 14, 'XTick', 1820:20:2020, 'YTick', -80:20:80);

    % Labels
    if seg == 1
        ylabel('Change in Variance (%)', 'FontSize', 16);
    else
        set(gca, 'YTickLabel', []);
    end
    if seg == 2
        xlabel('Year CE', 'FontSize', 16);
    end
    box on;
end

%% Legend 
axes(ax(end));
hold on;
h1 = plot(nan, nan, '-', 'Color', [0, 0, 0], 'LineWidth', 2);
h2 = plot(nan, nan, '-', 'Color', [0, 0, 1], 'LineWidth', 2);
h3 = plot(nan, nan, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);

legend([h1, h2, h3], ...
    {'Coral \delta^{18}O (Previously Published)', ...
     'Coral \delta^{18}O (This Study)', ...
     'Niño3.4 SST'}, ...
    'Location', 'northeast', 'FontSize', 14);
