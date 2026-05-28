% ============================================================
% Figure 4 release script
% Panels: Figure 4 statistics and wind-stress diagnostics
% Purpose: Compute forecast-phase SSH RMSE against SWOT and relate RMSE to wind forcing.
%
% Notes for release:
% - Run this script directly after confirming that the local data paths are available.
% - Local absolute paths are intentionally preserved from the submitted analysis version.
% - Keep figure settings unchanged unless you are intentionally updating the submitted figure.
% ============================================================

clearvars; close all; clc

%% ============================================================
% RMSE time series against SWOT
% 08-15 May 2023
% Simulations:
%   v6_fr_07_05  = GEN
%   v7_fr_07_05  = SWOT DA
%   v8_fr_07_05  = In-situ DA
%   v9_fr_07_05  = Joint DA
% X-axis uses actual SWOT observation time
%% ============================================================

addpath cmocean
addpath TEOS_10
addpath TEOS_10/library/

RT = load("rt_colormaps.mat");

%% ============================================================
% Paths
%% ============================================================
WMOP_path_v6 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v6_fr_07_05/'; % GEN
WMOP_path_v7 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/'; % SWOT DA
WMOP_path_v8 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8_fr_07_05/'; % In-situ DA
WMOP_path_v9 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'; % Joint DA

SWOT_path    = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';

WMOP_paths = {WMOP_path_v6, WMOP_path_v7, WMOP_path_v8, WMOP_path_v9};
labels = {'GEN','SWOT DA','In-situ DA','Joint DA'};
colors = [0.0660    0.4430    0.7450
          0.2310 0.6660 0.1960
          0.5210 0.0860 0.8190
          0.9608 0.4667 0.1608];

idx_joint = 4;   % Joint DA column in RMSE

dir_SWOT = dir(fullfile(SWOT_path,'*.nc'));

nSim = numel(WMOP_paths);
dir_mod = cell(nSim,1);
for s = 1:nSim
    dir_mod{s} = dir(fullfile(WMOP_paths{s},'*.nc'));
end

%% ============================================================
% Settings
%% ============================================================
target_dates  = datetime(2023,5,8):days(1):datetime(2023,5,15);

axis_def_mask = [0.8 3.0 39.0 41.0];
lon_min = axis_def_mask(1);
lon_max = axis_def_mask(2);
lat_min = axis_def_mask(3);
lat_max = axis_def_mask(4);

lon_target = 2.0;
lat_target = 40.0;

%% ============================================================
% Precompute model file day metadata
%% ============================================================
file_day = cell(nSim,1);

for s = 1:nSim
    this_dir = dir_mod{s};
    this_day = nan(length(this_dir),1);

    for k = 1:length(this_dir)
        try
            tt = ncread(fullfile(WMOP_paths{s}, this_dir(k).name), 'ocean_time');
            tt = datenum(seconds(tt) + datenum(1968,5,23));
            this_day(k) = floor(tt(1));
        catch
            this_day(k) = NaN;
        end
    end

    file_day{s} = this_day;
end

%% ============================================================
% Precompute SWOT representative times
% representative time = closest point in swath to (2E, 40N)
%% ============================================================
swot_closest_time = nan(length(dir_SWOT),1);
swot_day          = nan(length(dir_SWOT),1);

for k = 1:length(dir_SWOT)
    fname = fullfile(SWOT_path, dir_SWOT(k).name);
    try
        swot_time = ncread(fname,'time');
        swot_lon  = ncread(fname,'longitude');
        swot_lat  = ncread(fname,'latitude');

        TIME = datenum(seconds(swot_time) + datenum(2000,1,1));

        dist = sqrt((swot_lon - lon_target).^2 + (swot_lat - lat_target).^2);
        [min_dist_per_time, ~] = min(dist, [], 1);
        [~, best_time_idx] = min(min_dist_per_time);

        swot_closest_time(k) = TIME(best_time_idx);
        swot_day(k)          = floor(TIME(best_time_idx));
    catch
        swot_closest_time(k) = NaN;
        swot_day(k)          = NaN;
    end
end

%% ============================================================
% Containers
%% ============================================================
RMSE = nan(length(target_dates), nSim);
model_times_used = nan(length(target_dates), nSim);
swot_times_used  = nan(length(target_dates), 1);

%% ============================================================
% Loop over target days
%% ============================================================
for i = 1:length(target_dates)

    target_dnum = datenum(target_dates(i));

    %% --------------------------------------------------------
    % Find SWOT file for this day first
    %% --------------------------------------------------------
    idx_swot_day = find(swot_day == floor(target_dnum));

    if isempty(idx_swot_day)
        valid_swot = find(isfinite(swot_closest_time));
        [~, jj] = min(abs(swot_closest_time(valid_swot) - target_dnum));
        idx_swot = valid_swot(jj);
    else
        [~, jj] = min(abs(swot_closest_time(idx_swot_day) - target_dnum));
        idx_swot = idx_swot_day(jj);
    end

    [~, data_swot] = read_nc_file_struct(fullfile(SWOT_path, dir_SWOT(idx_swot).name));
    data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));

    t_swot = swot_closest_time(idx_swot);
    swot_times_used(i) = t_swot;

    mask_swot = data_swot.longitude >= lon_min & data_swot.longitude <= lon_max & ...
                data_swot.latitude  >= lat_min & data_swot.latitude  <= lat_max;

    ssha_masked   = data_swot.ssha_filtered;
    swot_lon_mask = data_swot.longitude;
    swot_lat_mask = data_swot.latitude;

    ssha_masked(~mask_swot)   = NaN;
    swot_lon_mask(~mask_swot) = NaN;
    swot_lat_mask(~mask_swot) = NaN;

    swot_ssh = (ssha_masked - nanmean(ssha_masked(:))) * 100;

    %% --------------------------------------------------------
    % Loop over simulations
    %% --------------------------------------------------------
    for s = 1:nSim

        this_file_day = file_day{s};
        this_dir      = dir_mod{s};
        this_path     = WMOP_paths{s};

        idx_day = find(this_file_day == floor(target_dnum));

        if isempty(idx_day)
            valid_idx = find(isfinite(this_file_day));
            [~, jj] = min(abs(this_file_day(valid_idx) - floor(target_dnum)));
            idx_mod = valid_idx(jj);
        else
            idx_mod = idx_day(1);
        end

        [~, data_mod] = read_nc_file_struct(fullfile(this_path, this_dir(idx_mod).name));
        data_mod.time = datenum(seconds(data_mod.ocean_time) + datenum(1968,5,23));

        [~, it_mod] = min(abs(data_mod.time - t_swot));
        model_times_used(i,s) = data_mod.time(it_mod);

        zeta = data_mod.zeta(:,:,it_mod)';
        [Lon, Lat] = meshgrid(data_mod.lon, data_mod.lat);

        mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;
        zeta_masked = zeta;
        zeta_masked(~mask) = NaN;
        zeta_anom = zeta - nanmean(zeta_masked(:));

        model_to_swot = interp2(data_mod.lon, data_mod.lat, zeta_anom*100, ...
                                swot_lon_mask(36:end,:), swot_lat_mask(36:end,:));

        A = swot_ssh(36:end,:);
        B = model_to_swot;
        RMSE(i,s) = sqrt(nanmean((A(:) - B(:)).^2));
    end
end

%% ============================================================
% Sort by actual SWOT observation time
%% ============================================================
[swot_times_used, isort] = sort(swot_times_used);
RMSE = RMSE(isort,:);
model_times_used = model_times_used(isort,:);

%% ============================================================
% ATM data
%% ============================================================
ATM_path_data  = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Atmospheric_data/';
[info_ATM, data_ATM] = read_nc_file_struct([ATM_path_data 'WMOP_hourly_wind_stress_heat_freshwater_flux_WMOP_REANALYSIS_FaSt_SWOT_v1.nc']);
data_ATM.time = datenum(seconds(data_ATM.ocean_time) + datenum(1968,05,23));

tau   = sqrt(data_ATM.sustr.^2 + data_ATM.svstr.^2);
tau_u = data_ATM.sustr;
tau_v = data_ATM.svstr;

lon = data_ATM.lon;
lat = data_ATM.lat;
time_atm = data_ATM.time;

x_target = 1.6;
y_target = 39.8;

dy = 10 / 111;
dx = 10 / (111 * cosd(y_target));

tau_u_point = nan(length(time_atm),1);
tau_v_point = nan(length(time_atm),1);

for t = 1:length(time_atm)
    Ufield = squeeze(data_ATM.sustr(:,:,t));
    Vfield = squeeze(data_ATM.svstr(:,:,t));

    ind_lon = lon >= (x_target - dx) & lon <= (x_target + dx);
    ind_lat = lat >= (y_target - dy) & lat <= (y_target + dy);

    tau_u_point(t) = mean(Ufield(ind_lon, ind_lat), 'all', 'omitnan');
    tau_v_point(t) = mean(Vfield(ind_lon, ind_lat), 'all', 'omitnan');
end

tau_point = sqrt(tau_u_point.^2 + tau_v_point.^2);

%% ============================================================
% Plot time series using actual SWOT times on x-axis
%% ============================================================
figure('Color','w','Position',[200 200 900 400]); hold on

t_swot = datetime(swot_times_used,'ConvertFrom','datenum');
t_atm  = datetime(time_atm,'ConvertFrom','datenum');

% yyaxis left
p = gobjects(1,nSim);

for s = 1:nSim
    p(s) = plot(t_swot, RMSE(:,s), '-o', ...
        'LineWidth', 2.2, ...
        'Color', colors(s,:), ...
        'MarkerFaceColor', colors(s,:), ...
        'MarkerSize', 7);
end

ylabel('RMSE (cm)')
xlim([datetime(datenum(2023,5,7),'ConvertFrom','datenum') ...
      datetime(swot_times_used(end),'ConvertFrom','datenum')])

ax = gca;
ax.YColor = 'k';
text(0.02, 0.92, '(g)', 'Units','normalized', 'FontWeight','bold', 'FontSize',14)

yyaxis right
p_tau = plot(t_atm, tau_point, '-','color',[.6 .6 .6], 'LineWidth', 1.5);
ylabel('|\tau| (N m^{-2})')

ax = gca;
ax.GridAlpha = 0.15;
ax.YColor = 'k';

set(gca,'FontSize',14,'LineWidth',1.3,'TickDir','out')
box on
grid on
xtickformat('dd-MMM')
xlim([datetime(datenum(2023,5,7),'ConvertFrom','datenum') ...
      datetime(swot_times_used(end),'ConvertFrom','datenum')])

legend([p, p_tau], [labels, {'|\tau|'}], ...
    'Location','northeast', 'Box','off')

legend([p], [labels], ...
    'Location','northeast', 'Box','off')

print(gcf, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_SWOT_stats', '-dpng', '-r300');

%% ============================================================
% RMSE summary table for panel (g)
%% ============================================================

% Simulation names
sim_names = {'GEN'; 'SWOT DA'; 'In-situ DA'; 'Joint DA'};

% If RMSE is already in cm, keep this:
RMSE_cm = RMSE;

% If RMSE was calculated in m, use this instead:
% RMSE_cm = RMSE * 100;

% Mean and standard deviation across SWOT dates
RMSE_Mean = mean(RMSE_cm, 1, 'omitnan')';
RMSE_Std  = std(RMSE_cm, 0, 1, 'omitnan')';

% Optional: improvement relative to GEN mean RMSE
RMSE_Improvement_vs_GEN = (RMSE_Mean(1) - RMSE_Mean) ./ RMSE_Mean(1) * 100;

% Create table
summary_RMSE_table = table( ...
    sim_names, ...
    RMSE_Mean, ...
    RMSE_Std, ...
    RMSE_Improvement_vs_GEN, ...
    'VariableNames', {'Simulation', 'SSH_RMSE_Mean_cm', 'SSH_RMSE_Std_cm', ...
                      'Improvement_vs_GEN_percent'});

disp(summary_RMSE_table)

% Optional: save table
writetable(summary_RMSE_table, 'summary_SSH_RMSE_table.csv');

return

% % ============================================================
% % Figure 4 release script
% % Panels: Figure 4g
% % Purpose: Forecast-period SSH RMSE time series against SWOT, with wind-stress comparison.
% %
% % Notes for release:
% % - Run this script for panel (g).
% % - It uses forecast simulations and SWOT overpasses from 8-15 May.
% % - This is separate from the Figure 4 map and ADCP scripts.
% % ============================================================
%
%
% clear; close all; clc
%
% %% ============================================================
% % RMSE time series against SWOT
% % 08-15 May 2023
% % Simulations:
% %   v7_fr_07_05  = SWOT DA
% %   v8_fr_07_05  = In-situ DA
% %   v9_fr_07_05  = Joint DA
% % X-axis uses actual SWOT observation time
% %% ============================================================
%
% addpath cmocean
% addpath TEOS_10
% addpath TEOS_10/library/
%
% RT = load("rt_colormaps.mat");
%
% %% ============================================================
% % Paths
% %% ============================================================
% WMOP_path_v7 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/'; % SWOT DA
% WMOP_path_v8 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8_fr_07_05/'; % In-situ DA
% WMOP_path_v9 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'; % Joint DA
%
% SWOT_path    = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';
%
% dir_v7 = dir(fullfile(WMOP_path_v7,'*.nc'));
% dir_v8 = dir(fullfile(WMOP_path_v8,'*.nc'));
% dir_v9 = dir(fullfile(WMOP_path_v9,'*.nc'));
% dir_SWOT = dir(fullfile(SWOT_path,'*.nc'));
%
% %% ============================================================
% % Settings
% %% ============================================================
% target_dates  = datetime(2023,5,8):days(1):datetime(2023,5,15);
%
% axis_def_mask = [0.8 3.0 39.0 41.0];
% lon_min = axis_def_mask(1);
% lon_max = axis_def_mask(2);
% lat_min = axis_def_mask(3);
% lat_max = axis_def_mask(4);
%
% lon_target = 2.0;
% lat_target = 40.0;
%
% colors =[ 0.2310    0.6660    0.1960
%           0.5210    0.0860    0.8190
%           0.9608    0.4667    0.1608 ];
%
% labels = {'SWOT DA','In-situ DA','Joint DA'};
%
% %% ============================================================
% % Precompute model file day metadata
% %% ============================================================
% file_day_v7 = nan(length(dir_v7),1);
% file_day_v8 = nan(length(dir_v8),1);
% file_day_v9 = nan(length(dir_v9),1);
%
% for k = 1:length(dir_v7)
%     try
%         tt = ncread(fullfile(WMOP_path_v7,dir_v7(k).name),'ocean_time');
%         tt = datenum(seconds(tt) + datenum(1968,5,23));
%         file_day_v7(k) = floor(tt(1));
%     catch
%         file_day_v7(k) = NaN;
%     end
% end
%
% for k = 1:length(dir_v8)
%     try
%         tt = ncread(fullfile(WMOP_path_v8,dir_v8(k).name),'ocean_time');
%         tt = datenum(seconds(tt) + datenum(1968,5,23));
%         file_day_v8(k) = floor(tt(1));
%     catch
%         file_day_v8(k) = NaN;
%     end
% end
%
% for k = 1:length(dir_v9)
%     try
%         tt = ncread(fullfile(WMOP_path_v9,dir_v9(k).name),'ocean_time');
%         tt = datenum(seconds(tt) + datenum(1968,5,23));
%         file_day_v9(k) = floor(tt(1));
%     catch
%         file_day_v9(k) = NaN;
%     end
% end
%
% %% ============================================================
% % Precompute SWOT representative times
% % representative time = closest point in swath to (2E, 40N)
% %% ============================================================
% swot_closest_time = nan(length(dir_SWOT),1);
% swot_day          = nan(length(dir_SWOT),1);
%
% for k = 1:length(dir_SWOT)
%     fname = fullfile(SWOT_path, dir_SWOT(k).name);
%     try
%         swot_time = ncread(fname,'time');
%         swot_lon  = ncread(fname,'longitude');
%         swot_lat  = ncread(fname,'latitude');
%
%         TIME = datenum(seconds(swot_time) + datenum(2000,1,1));
%
%         dist = sqrt((swot_lon - lon_target).^2 + (swot_lat - lat_target).^2);
%         [min_dist_per_time, ~] = min(dist, [], 1);
%         [~, best_time_idx] = min(min_dist_per_time);
%
%         swot_closest_time(k) = TIME(best_time_idx);
%         swot_day(k)          = floor(TIME(best_time_idx));
%     catch
%         swot_closest_time(k) = NaN;
%         swot_day(k)          = NaN;
%     end
% end
%
% %% ============================================================
% % Containers
% %% ============================================================
% RMSE = nan(length(target_dates),3);
% model_times_used = nan(length(target_dates),3);
% swot_times_used  = nan(length(target_dates),1);
%
% %% ============================================================
% % Loop over target days
% %% ============================================================
% for i = 1:length(target_dates)
%
%     target_dnum = datenum(target_dates(i));
%
%     %% --------------------------------------------------------
%     % Find SWOT file for this day first
%     %% --------------------------------------------------------
%     idx_swot_day = find(swot_day == floor(target_dnum));
%
%     if isempty(idx_swot_day)
%         valid_swot = find(isfinite(swot_closest_time));
%         [~, jj] = min(abs(swot_closest_time(valid_swot) - target_dnum));
%         idx_swot = valid_swot(jj);
%     else
%         [~, jj] = min(abs(swot_closest_time(idx_swot_day) - target_dnum));
%         idx_swot = idx_swot_day(jj);
%     end
%
%     [~, data_swot] = read_nc_file_struct(fullfile(SWOT_path, dir_SWOT(idx_swot).name));
%     data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));
%
%     t_swot = swot_closest_time(idx_swot);
%     swot_times_used(i) = t_swot;
%
%     mask_swot = data_swot.longitude >= lon_min & data_swot.longitude <= lon_max & ...
%                 data_swot.latitude  >= lat_min & data_swot.latitude  <= lat_max;
%
%     ssha_masked   = data_swot.ssha_filtered;
%     swot_lon_mask = data_swot.longitude;
%     swot_lat_mask = data_swot.latitude;
%
%     ssha_masked(~mask_swot)   = NaN;
%     swot_lon_mask(~mask_swot) = NaN;
%     swot_lat_mask(~mask_swot) = NaN;
%
%     swot_ssh = (ssha_masked - nanmean(ssha_masked(:))) * 100;
%
%     %% ---------------- SWOT DA ----------------
%     idx_day = find(file_day_v7 == floor(target_dnum));
%     if isempty(idx_day)
%         valid_idx = find(isfinite(file_day_v7));
%         [~, jj] = min(abs(file_day_v7(valid_idx) - floor(target_dnum)));
%         idx_mod = valid_idx(jj);
%     else
%         idx_mod = idx_day(1);
%     end
%
%     [~, data_mod] = read_nc_file_struct(fullfile(WMOP_path_v7, dir_v7(idx_mod).name));
%     data_mod.time = datenum(seconds(data_mod.ocean_time) + datenum(1968,5,23));
%
%     [~, it_mod] = min(abs(data_mod.time - t_swot));
%     model_times_used(i,1) = data_mod.time(it_mod);
%
%     zeta = data_mod.zeta(:,:,it_mod)';
%     [Lon, Lat] = meshgrid(data_mod.lon, data_mod.lat);
%     mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;
%     zeta_masked = zeta;
%     zeta_masked(~mask) = NaN;
%     zeta_anom = zeta - nanmean(zeta_masked(:));
%
%     model_to_swot = interp2(data_mod.lon, data_mod.lat, zeta_anom*100, ...
%                             swot_lon_mask(36:end,:), swot_lat_mask(36:end,:));
%
%     A = swot_ssh(36:end,:);
%     B = model_to_swot;
%     RMSE(i,1) = sqrt(nanmean((A(:) - B(:)).^2));
%
%     %% ---------------- In-situ DA ----------------
%     idx_day = find(file_day_v8 == floor(target_dnum));
%     if isempty(idx_day)
%         valid_idx = find(isfinite(file_day_v8));
%         [~, jj] = min(abs(file_day_v8(valid_idx) - floor(target_dnum)));
%         idx_mod = valid_idx(jj);
%     else
%         idx_mod = idx_day(1);
%     end
%
%     [~, data_mod] = read_nc_file_struct(fullfile(WMOP_path_v8, dir_v8(idx_mod).name));
%     data_mod.time = datenum(seconds(data_mod.ocean_time) + datenum(1968,5,23));
%
%     [~, it_mod] = min(abs(data_mod.time - t_swot));
%     model_times_used(i,2) = data_mod.time(it_mod);
%
%     zeta = data_mod.zeta(:,:,it_mod)';
%     [Lon, Lat] = meshgrid(data_mod.lon, data_mod.lat);
%     mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;
%     zeta_masked = zeta;
%     zeta_masked(~mask) = NaN;
%     zeta_anom = zeta - nanmean(zeta_masked(:));
%
%     model_to_swot = interp2(data_mod.lon, data_mod.lat, zeta_anom*100, ...
%                             swot_lon_mask(36:end,:), swot_lat_mask(36:end,:));
%
%     A = swot_ssh(36:end,:);
%     B = model_to_swot;
%     RMSE(i,2) = sqrt(nanmean((A(:) - B(:)).^2));
%
%     %% ---------------- Joint DA ----------------
%     idx_day = find(file_day_v9 == floor(target_dnum));
%     if isempty(idx_day)
%         valid_idx = find(isfinite(file_day_v9));
%         [~, jj] = min(abs(file_day_v9(valid_idx) - floor(target_dnum)));
%         idx_mod = valid_idx(jj);
%     else
%         idx_mod = idx_day(1);
%     end
%
%     [~, data_mod] = read_nc_file_struct(fullfile(WMOP_path_v9, dir_v9(idx_mod).name));
%     data_mod.time = datenum(seconds(data_mod.ocean_time) + datenum(1968,5,23));
%
%     [~, it_mod] = min(abs(data_mod.time - t_swot));
%     model_times_used(i,3) = data_mod.time(it_mod);
%
%     zeta = data_mod.zeta(:,:,it_mod)';
%     [Lon, Lat] = meshgrid(data_mod.lon, data_mod.lat);
%     mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;
%     zeta_masked = zeta;
%     zeta_masked(~mask) = NaN;
%     zeta_anom = zeta - nanmean(zeta_masked(:));
%
%     model_to_swot = interp2(data_mod.lon, data_mod.lat, zeta_anom*100, ...
%                             swot_lon_mask(36:end,:), swot_lat_mask(36:end,:));
%
%     A = swot_ssh(36:end,:);
%     B = model_to_swot;
%     RMSE(i,3) = sqrt(nanmean((A(:) - B(:)).^2));
% end
%
% %% ============================================================
% % Sort by actual SWOT observation time
% %% ============================================================
% [swot_times_used, isort] = sort(swot_times_used);
% RMSE = RMSE(isort,:);
% model_times_used = model_times_used(isort,:);
%
% %% ============================================================
% % Plot time series using actual SWOT times on x-axis
% %% ============================================================
%
%
%
% %% ATM data
% ATM_path_data  = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Atmospheric_data/';
% [info_ATM, data_ATM] = read_nc_file_struct([ATM_path_data 'WMOP_hourly_wind_stress_heat_freshwater_flux_WMOP_REANALYSIS_FaSt_SWOT_v1.nc']);
% data_ATM.time = datenum(seconds(data_ATM.ocean_time) + datenum(1968,05,23));
% %% Precompute wind quantities
% tau   = sqrt(data_ATM.sustr.^2 + data_ATM.svstr.^2);
% tau_u = data_ATM.sustr;
% tau_v = data_ATM.svstr;
% wind_dir = atan2d(tau_v, tau_u); % wind direction in degrees (-180 to 180)
% wind_time = data_ATM.time;
%
% %% Extract lon, lat, time from ATM data
% lon = data_ATM.lon;   % assuming size (nx)
% lat = data_ATM.lat;   % assuming size (ny)
% time_atm = data_ATM.time;  % already in datenum
%
% lon = data_ATM.lon;   % size (nx)
% lat = data_ATM.lat;   % size (ny)
% [LON, LAT] = meshgrid(lon, lat);   % 2D mesh for quiver
%
%
%
% % Target point
% x_target = 1.6;
% y_target = 39.8;
%
% % x_target = 5;
% % y_target = 40.5;
%
% % Convert 10 km to degrees
% dy = 10 / 111;
% dx = 10 / (111 * cosd(y_target));
%
% % Preallocate
% tau_u_point = nan(length(time_atm),1);
% tau_v_point = nan(length(time_atm),1);
%
% for t = 1:length(time_atm)
%     % Extract 2D fields at time t
%     Ufield = squeeze(data_ATM.sustr(:,:,t));
%     Vfield = squeeze(data_ATM.svstr(:,:,t));
%     shfluxfield = squeeze(data_ATM.shflux(:,:,t));
%     EmPfield    = squeeze(data_ATM.EminusP(:,:,t));
%
%     % Find grid indices within ±dx, ±dy around target
%     ind_lon = lon >= (x_target - dx) & lon <= (x_target + dx);
%     ind_lat = lat >= (y_target - dy) & lat <= (y_target + dy);
%
%     % Average over box
%     tau_u_point(t) = mean(Ufield(ind_lon, ind_lat), 'all', 'omitnan');
%     tau_v_point(t) = mean(Vfield(ind_lon, ind_lat), 'all', 'omitnan');
%     % shflux(t)      = mean(shfluxfield(ind_lon, ind_lat), 'all', 'omitnan');
%     % EmP(t)         = mean(EmPfield(ind_lon, ind_lat), 'all', 'omitnan');
% end
%
%
% % return
% % Magnitude and direction
% tau_point = sqrt(tau_u_point.^2 + tau_v_point.^2);
%
%
%
% figure('Color','w','Position',[200 200 850 380]); hold on
%
% t_swot = datetime(swot_times_used,'ConvertFrom','datenum');
% t_atm  = datetime(time_atm,'ConvertFrom','datenum');   % only if time_atm is datenum
%
% yyaxis left
% p1 = plot(t_swot, RMSE(:,1), '-o', ...
%     'LineWidth', 2.2, 'Color', colors(1,:), ...
%     'MarkerFaceColor', colors(1,:), 'MarkerSize', 7);
%
% p2 = plot(t_swot, RMSE(:,2), '-o', ...
%     'LineWidth', 2.2, 'Color', colors(2,:), ...
%     'MarkerFaceColor', colors(2,:), 'MarkerSize', 7);
%
% p3 = plot(t_swot, RMSE(:,3), '-o', ...
%     'LineWidth', 2.2, 'Color', colors(3,:), ...
%     'MarkerFaceColor', colors(3,:), 'MarkerSize', 7);
%
% ylabel('RMSE (cm)')
% xlim([datetime(datenum(2023,5,8),'ConvertFrom','datenum') ...
%       datetime(swot_times_used(end),'ConvertFrom','datenum')])
%
% ax = gca;
% ax.YColor = 'k';
% text(0.02, 0.92, '(g)', 'Units','normalized', 'FontWeight','bold', 'FontSize',14)
%
%
% yyaxis right
% p4 = plot(t_atm, tau_point, 'k-', 'LineWidth', 1.5);
% ylabel('|\tau| (N m^{-2})')
%
% set(gca,'FontSize',14,'LineWidth',1.3,'TickDir','out')
% box on
% grid on
% ax = gca;
% ax.GridAlpha = 0.15;
%
% xtickformat('dd-MMM')
% % xlim([datetime(2023,5,8) datetime(2023,5,15)])
% xlim([datetime(datenum(2023,5,7),'ConvertFrom','datenum') ...
%       datetime(swot_times_used(end),'ConvertFrom','datenum')])
%
% ax = gca;
% ax.YColor = 'k';
%
% legend([p1,p2,p3,p4], [labels, {'|\tau|'}], ...
%     'Location','northeast', 'Box','off')
%
% print(gcf, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_SWOT_stats_wind', '-dpng', '-r300');
%
%
% return
% figure
% p4 = plot(t_atm, tau_point, 'k-', 'LineWidth', 1.5);
% xlim([datetime(datenum(2023,4,21),'ConvertFrom','datenum') ...
%       datetime(datenum(2023,5,6),'ConvertFrom','datenum')])
%
%
%
% %% Joint DA RMSE vs maximum wind stress over previous 6 h
% % Assumes variables already exist:
% % t_swot, t_atm, RMSE, tau_point
%
% rmse_joint = RMSE(:,3);
%
% % Convert time to datenum if needed
% if isdatetime(t_swot)
%     x_swot = datenum(t_swot);
% else
%     x_swot = t_swot;
% end
%
% if isdatetime(t_atm)
%     x_atm = datenum(t_atm);
% else
%     x_atm = t_atm;
% end
%
% % Compute max tau over [t-6h, t] for each RMSE time
% tau_6h_max = nan(size(x_swot));
%
% for i = 1:length(x_swot)
%     ind = x_atm >= (x_swot(i) - 6/24) & x_atm <= x_swot(i);
%     if any(ind)
%         tau_6h_max(i) = max(tau_point(ind), [], 'omitnan');
%     end
% end
%
% % Optional: restrict to prediction period only
% use_prediction_only = true;
%
% if use_prediction_only
%     if isdatetime(t_swot)
%         t1 = datetime(2023,5,7);
%         t2 = datetime(2023,5,15,23,59,59);
%         mask = t_swot >= t1 & t_swot <= t2;
%     else
%         t1 = datenum(2023,5,7);
%         t2 = datenum(2023,5,15,23,59,59);
%         mask = x_swot >= t1 & x_swot <= t2;
%     end
% else
%     mask = true(size(x_swot));
% end
%
% % Remove NaNs
% ok = mask & isfinite(rmse_joint) & isfinite(tau_6h_max);
%
% x = tau_6h_max(ok);
% y = rmse_joint(ok);
%
% % Pearson correlation
% [r, p] = corr(x(:), y(:), 'Type', 'Pearson');
%
% % Linear fit
% coef = polyfit(x, y, 1);
% xfit = linspace(min(x), max(x), 100);
% yfit = polyval(coef, xfit);
%
% % Plot
% figure('Color','w','Position',[200 150 620 500]);
% hold on; box on
%
% scatter(x, y, 90, 'filled', ...
%     'MarkerFaceColor', [0.85 0.33 0.10], ...
%     'MarkerEdgeColor', 'k', ...
%     'LineWidth', 0.8);
%
% plot(xfit, yfit, 'k-', 'LineWidth', 2.0);
%
% xlabel('Max wind stress over previous 6 h (N m^{-2})', 'FontSize', 16)
% ylabel('Joint DA RMSE (cm)', 'FontSize', 16)
%
% set(gca, 'FontSize', 15, ...
%          'LineWidth', 1.2, ...
%          'TickDir', 'out', ...
%          'Layer', 'top')
%
% grid on
%
% txt = sprintf('r = %.2f\np = %.3f\nn = %d', r, p, numel(x));
% text(0.05, 0.95, txt, 'Units', 'normalized', ...
%     'FontSize', 14, ...
%     'VerticalAlignment', 'top', ...
%     'BackgroundColor', 'w', ...
%     'Margin', 6);
%
% title('Joint DA RMSE vs max wind stress over previous 6 h', 'FontSize', 17)
%
% fprintf('Pearson r = %.3f, p = %.4f, n = %d\n', r, p, numel(x));
% fprintf('Linear fit: RMSE = %.3f * tau_{6h,max} + %.3f\n', coef(1), coef(2));
%
% % Optional: show paired values
% if isdatetime(t_swot)
%     T = table(t_swot(ok), x, y, ...
%         'VariableNames', {'Time','Tau_6h_max','RMSE_JointDA'});
% else
%     T = table(x_swot(ok), x, y, ...
%         'VariableNames', {'Time_datenum','Tau_6h_max','RMSE_JointDA'});
% end
% disp(T)
%
%
%
%
%
% %% RMSE vs preceding 6 h mean wind stress
% % Assumes variables already exist:
% % t_swot, t_atm, RMSE, tau_point
%
% rmse_joint = RMSE(:,3);
%
% % Convert time to datenum if needed
% if isdatetime(t_swot)
%     x_swot = datenum(t_swot);
% else
%     x_swot = t_swot;
% end
%
% if isdatetime(t_atm)
%     x_atm = datenum(t_atm);
% else
%     x_atm = t_atm;
% end
%
% % Compute mean tau over [t-6h, t] for each RMSE time
% tau_6h = nan(size(x_swot));
%
% for i = 1:length(x_swot)
%     ind = x_atm >= (x_swot(i) - 3/24) & x_atm <= x_swot(i);
%     tau_6h(i) = mean(tau_point(ind), 'omitnan');
% end
%
% % Optional: restrict to prediction period only
% use_prediction_only = true;
%
% if use_prediction_only
%     if isdatetime(t_swot)
%         t1 = datetime(2023,5,7);
%         t2 = datetime(2023,5,15,23,59,59);
%         mask = t_swot >= t1 & t_swot <= t2;
%     else
%         t1 = datenum(2023,5,7);
%         t2 = datenum(2023,5,15,23,59,59);
%         mask = x_swot >= t1 & x_swot <= t2;
%     end
% else
%     mask = true(size(x_swot));
% end
%
% % Remove NaNs
% ok = mask & isfinite(rmse_joint) & isfinite(tau_6h);
%
% x = tau_6h(ok);
% y = rmse_joint(ok);
%
% % Pearson correlation
% [r, p] = corr(x(:), y(:), 'Type', 'Pearson');
%
% % Linear fit
% coef = polyfit(x, y, 1);
% xfit = linspace(min(x), max(x), 100);
% yfit = polyval(coef, xfit);
%
% % Plot
% figure('Color','w','Position',[200 150 620 500]);
% hold on; box on
%
% scatter(x, y, 90, 'filled', ...
%     'MarkerFaceColor', [0.85 0.33 0.10], ...
%     'MarkerEdgeColor', 'k', ...
%     'LineWidth', 0.8);
%
% plot(xfit, yfit, 'k-', 'LineWidth', 2.0);
%
% xlabel('Mean wind stress over previous 6 h (N m^{-2})', 'FontSize', 16)
% ylabel('Joint DA RMSE (cm)', 'FontSize', 16)
%
% set(gca, 'FontSize', 15, ...
%          'LineWidth', 1.2, ...
%          'TickDir', 'out', ...
%          'Layer', 'top')
%
% grid on
%
% txt = sprintf('r = %.2f\np = %.3f\nn = %d', r, p, numel(x));
% text(0.05, 0.95, txt, 'Units', 'normalized', ...
%     'FontSize', 14, ...
%     'VerticalAlignment', 'top', ...
%     'BackgroundColor', 'w', ...
%     'Margin', 6);
%
% title('Joint DA RMSE vs mean wind stress over previous 6 h', 'FontSize', 17)
%
% fprintf('Pearson r = %.3f, p = %.4f, n = %d\n', r, p, numel(x));
% fprintf('Linear fit: RMSE = %.3f * tau_6h + %.3f\n', coef(1), coef(2));
%
% % Optional: display paired values
% if isdatetime(t_swot)
%     T = table(t_swot(ok), x, y, ...
%         'VariableNames', {'Time','Tau_6h_mean','RMSE_JointDA'});
% else
%     T = table(x_swot(ok), x, y, ...
%         'VariableNames', {'Time_datenum','Tau_6h_mean','RMSE_JointDA'});
% end
% disp(T)
%
%
%
%
%
%
%
% %% Joint DA RMSE vs wind stress: Pearson r
% rmse_joint = RMSE(:,3);
%
% % Convert time to numeric if needed
% if isdatetime(t_swot)
%     x_swot = datenum(t_swot);
% else
%     x_swot = t_swot;
% end
%
% if isdatetime(t_atm)
%     x_atm = datenum(t_atm);
% else
%     x_atm = t_atm;
% end
%
% % Interpolate wind stress to the RMSE times
% tau_match = interp1(x_atm, tau_point, x_swot, 'linear');
%
% % Remove NaNs
% ok = isfinite(rmse_joint) & isfinite(tau_match);
%
% % Pearson correlation
% [r, p] = corr(tau_match(ok), rmse_joint(ok), 'Type', 'Pearson');
%
% fprintf('Joint DA vs tau at matched times: r = %.3f, p = %.4f, n = %d\n', ...
%         r, p, sum(ok));
%
% % Optional: show paired values
% T = table(t_swot(ok), rmse_joint(ok), tau_match(ok), ...
%     'VariableNames', {'Time','RMSE_JointDA','Tau_matched'});
% disp(T)
%
%
%
%
% %% Scatter plot with normalized x and y
% % x = wind stress at RMSE times
% % y = Joint DA RMSE
% % Variables assumed to exist:
% % t_swot, t_atm, RMSE, tau_point
%
% rmse_joint = RMSE(:,3);
%
% % Convert time to numeric if needed
% if isdatetime(t_swot)
%     x_swot = datenum(t_swot);
% else
%     x_swot = t_swot;
% end
%
% if isdatetime(t_atm)
%     x_atm = datenum(t_atm);
% else
%     x_atm = t_atm;
% end
%
% % -------- choose wind metric --------
% % Option 1: instantaneous tau at RMSE times
% tau_match = interp1(x_atm, tau_point, x_swot, 'linear');
%
% % % Option 2: preceding 24 h mean tau, more physical
% % tau_match = nan(size(x_swot));
% % for i = 1:length(x_swot)
% %     ind = x_atm >= (x_swot(i)-1) & x_atm <= x_swot(i);
% %     tau_match(i) = mean(tau_point(ind), 'omitnan');
% % end
%
% % -------- choose period --------
% use_prediction_only = true;
%
% if use_prediction_only
%     if isdatetime(t_swot)
%         t1 = datetime(2023,5,7);
%         t2 = datetime(2023,5,15,23,59,59);
%         mask = t_swot >= t1 & t_swot <= t2;
%     else
%         t1 = datenum(2023,5,7);
%         t2 = datenum(2023,5,15,23,59,59);
%         mask = x_swot >= t1 & x_swot <= t2;
%     end
% else
%     mask = true(size(rmse_joint));
% end
%
% % Remove NaNs
% ok = mask & isfinite(rmse_joint) & isfinite(tau_match);
%
% x = tau_match(ok);
% y = rmse_joint(ok);
%
% % -------- normalize both x and y using z-score --------
% x_norm = (x - mean(x,'omitnan')) ./ std(x,'omitnan');
% y_norm = (y - mean(y,'omitnan')) ./ std(y,'omitnan');
%
% % Pearson correlation
% [r, p] = corr(x_norm, y_norm, 'Type', 'Pearson');
%
% % Linear fit in normalized space
% coef = polyfit(x_norm, y_norm, 1);
% xfit = linspace(min(x_norm), max(x_norm), 100);
% yfit = polyval(coef, xfit);
%
% % -------- plot --------
% figure('Color','w','Position',[200 150 620 500]);
% hold on; box on
%
% scatter(x_norm, y_norm, 85, 'filled', ...
%     'MarkerFaceColor', [0.85 0.33 0.10], ...
%     'MarkerEdgeColor', 'k', ...
%     'LineWidth', 0.8);
%
% plot(xfit, yfit, 'k-', 'LineWidth', 2.0);
%
% xlabel('Normalized wind stress', 'FontSize', 16)
% ylabel('Normalized Joint DA RMSE', 'FontSize', 16)
%
% set(gca, 'FontSize', 15, ...
%          'LineWidth', 1.2, ...
%          'TickDir', 'out', ...
%          'Layer', 'top')
%
% grid on
%
% txt = sprintf('r = %.2f\\np = %.3f\\nn = %d', r, p, numel(x_norm));
% text(0.05, 0.95, txt, 'Units', 'normalized', ...
%     'FontSize', 14, ...
%     'VerticalAlignment', 'top', ...
%     'BackgroundColor', 'w', ...
%     'Margin', 6);
%
% title('Normalized Joint DA RMSE vs wind stress', 'FontSize', 17)
%
% fprintf('Normalized Pearson r = %.3f, p = %.4f, n = %d\n', r, p, numel(x_norm));
% fprintf('Fit in normalized space: y = %.3f x + %.3f\n', coef(1), coef(2));
%
%
%
%
%
%
%
% %% Scatter plot: Joint DA RMSE vs wind stress
% % Variables assumed to exist:
% % t_swot, t_atm, RMSE, tau_point
%
% rmse_joint = RMSE(:,3);
%
% % Convert time to numeric if needed
% if isdatetime(t_swot)
%     x_swot = datenum(t_swot);
% else
%     x_swot = t_swot;
% end
%
% if isdatetime(t_atm)
%     x_atm = datenum(t_atm);
% else
%     x_atm = t_atm;
% end
%
% % Interpolate wind stress to RMSE times
% tau_match = interp1(x_atm, tau_point, x_swot, 'linear');
%
% % Optional: restrict to prediction period only
% use_prediction_only = true;
%
% if use_prediction_only
%     if isdatetime(t_swot)
%         t1 = datetime(2023,5,7);
%         t2 = datetime(2023,5,15,23,59,59);
%         mask = t_swot >= t1 & t_swot <= t2;
%     else
%         t1 = datenum(2023,5,7);
%         t2 = datenum(2023,5,15,23,59,59);
%         mask = x_swot >= t1 & x_swot <= t2;
%     end
% else
%     mask = true(size(rmse_joint));
% end
%
% % Remove NaNs
% ok = mask & isfinite(rmse_joint) & isfinite(tau_match);
%
% x = tau_match(ok);
% y = rmse_joint(ok);
%
% % Pearson correlation
% [r, p] = corr(x, y, 'Type', 'Pearson');
%
% % Linear fit
% coef = polyfit(x, y, 1);
% xfit = linspace(min(x), max(x), 100);
% yfit = polyval(coef, xfit);
%
% % Plot
% figure('Color','w','Position',[200 150 620 500]);
% hold on; box on
%
% scatter(x, y, 80, 'filled', ...
%     'MarkerFaceColor', [0.85 0.33 0.10], ...
%     'MarkerEdgeColor', 'k', ...
%     'LineWidth', 0.8);
%
% plot(xfit, yfit, 'k-', 'LineWidth', 2.0);
%
% xlabel('Wind stress magnitude |\tau| (N m^{-2})', 'FontSize', 16)
% ylabel('Joint DA RMSE (cm)', 'FontSize', 16)
%
% set(gca, 'FontSize', 15, 'LineWidth', 1.2, ...
%     'TickDir', 'out', 'Layer', 'top')
%
% grid on
%
% % Annotation
% txt = sprintf('r = %.2f\np = %.3f\nn = %d', r, p, numel(x));
% text(0.05, 0.95, txt, 'Units', 'normalized', ...
%     'FontSize', 14, 'VerticalAlignment', 'top', ...
%     'BackgroundColor', 'w', 'Margin', 6);
%
% title('Joint DA RMSE vs wind stress', 'FontSize', 17)
%
% fprintf('Pearson r = %.3f, p = %.4f, n = %d\n', r, p, numel(x));
% fprintf('Linear fit: RMSE = %.3f * tau + %.3f\n', coef(1), coef(2));
