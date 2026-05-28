% ============================================================
% Supporting Information release script
% Purpose: Plot daily temperature and salinity evolution at approximately 300 m from SWOT DA.
% Outputs: Supporting figures for cold-water and salinity evolution.
%
% Notes for release:
% - Run this script directly after confirming that the local data paths are available.
% - Local absolute paths are intentionally preserved from the submitted analysis version.
% - Keep figure settings unchanged unless you are intentionally updating the submitted figure.
% ============================================================

clearvars; close all; clc

%% ============================================================
% Daily evolution at ~300 m from SWOT DA
% Figure 1: Temperature + velocity vectors
% Figure 2: Salinity + velocity vectors
% ============================================================

addpath cmocean
RT = load("rt_colormaps.mat");

C = load('coastline_full_westmed_nolakes.mat');

ncst = C.ncst;

%% ============================================================
% Paths
% ============================================================
WMOP_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7/'; % SWOT DA
out_path  = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/';

dir_WMOP = dir(fullfile(WMOP_path, '*.nc'));

%% ============================================================
% Settings
% ============================================================
target_dates = datetime(2023,4,21):days(1):datetime(2023,5,6);
nDays = numel(target_dates);

hour_id = 6;                 % target hourly index
target_depth = 300;          % m
axis_def = [0.5 3 39 41];    % [lon_min lon_max lat_min lat_max]

lon_min = axis_def(1);
lon_max = axis_def(2);
lat_min = axis_def(3);
lat_max = axis_def(4);

% quiver thinning
skip_x = 6;
skip_y = 6;
quiv_scale = 2.5;

% subplot layout
ncol = 4;
nrow = 4;

%% ============================================================
% Read file times first and find matching daily files
% ============================================================
file_dates = NaT(length(dir_WMOP),1);

for ii = 1:length(dir_WMOP)
    fname = fullfile(WMOP_path, dir_WMOP(ii).name);

    try
        ocean_time = ncread(fname, 'ocean_time');
        t = datenum(seconds(ocean_time) + datenum(1968,5,23));
        file_dates(ii) = datetime(t(1), 'ConvertFrom', 'datenum');
    catch
        file_dates(ii) = NaT;
    end
end

file_days = dateshift(file_dates, 'start', 'day');

match_idx = nan(nDays,1);
for ii = 1:nDays
    kk = find(file_days == target_dates(ii), 1, 'first');
    if ~isempty(kk)
        match_idx(ii) = kk;
    end
end

if any(isnan(match_idx))
    warning('Some dates were not found in the WMOP files.')
    disp(table(target_dates(:), match_idx, 'VariableNames', {'TargetDate','FileIndex'}))
end

%% ============================================================
% Pre-read one file to get grid and depth index
% ============================================================
first_valid = find(~isnan(match_idx), 1, 'first');
if isempty(first_valid)
    error('No matching files found for the requested dates.')
end

sample_file = fullfile(WMOP_path, dir_WMOP(match_idx(first_valid)).name);

lon   = double(ncread(sample_file, 'lon'));
lat   = double(ncread(sample_file, 'lat'));
depth = double(ncread(sample_file, 'depth'));

[~, iz] = min(abs(depth - target_depth));
depth_used = depth(iz);

[Lon, Lat] = meshgrid(lon, lat);
mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;

ix = find(lon >= lon_min & lon <= lon_max);
iy = find(lat >= lat_min & lat <= lat_max);

lon_sub = lon(ix);
lat_sub = lat(iy);

%% ============================================================
% Read all daily fields first so color limits are consistent
% ============================================================
T_all = cell(nDays,1);
S_all = cell(nDays,1);
U_all = cell(nDays,1);
V_all = cell(nDays,1);
time_str = strings(nDays,1);

for ii = 1:nDays
    if isnan(match_idx(ii))
        continue
    end

    fname = fullfile(WMOP_path, dir_WMOP(match_idx(ii)).name);

    ocean_time = ncread(fname, 'ocean_time');
    t = datenum(seconds(ocean_time) + datenum(1968,5,23));
    model_time = datetime(t(hour_id), 'ConvertFrom', 'datenum');
    time_str(ii) = string(datestr(t(hour_id), 'dd-mmm'));

    % NOTE:
    % assumed variable dimensions are (lon, lat, depth, time)
    T = squeeze(ncread(fname, 'temp', [1 1 iz hour_id], [Inf Inf 1 1]))';
    S = squeeze(ncread(fname, 'salt', [1 1 iz hour_id], [Inf Inf 1 1]))';
    U = squeeze(ncread(fname, 'u',    [1 1 iz hour_id], [Inf Inf 1 1]))';
    V = squeeze(ncread(fname, 'v',    [1 1 iz hour_id], [Inf Inf 1 1]))';

    T(~mask) = NaN;
    S(~mask) = NaN;
    U(~mask) = NaN;
    V(~mask) = NaN;

    T_all{ii} = T(iy, ix);
    S_all{ii} = S(iy, ix);
    U_all{ii} = U(iy, ix);
    V_all{ii} = V(iy, ix);
end

% color limits from all valid days
Tmin = inf; Tmax = -inf;
Smin = inf; Smax = -inf;

for ii = 1:nDays
    if ~isempty(T_all{ii})
        Tmin = min(Tmin, min(T_all{ii}(:), [], 'omitnan'));
        Tmax = max(Tmax, max(T_all{ii}(:), [], 'omitnan'));
    end
    if ~isempty(S_all{ii})
        Smin = min(Smin, min(S_all{ii}(:), [], 'omitnan'));
        Smax = max(Smax, max(S_all{ii}(:), [], 'omitnan'));
    end
end

%% ============================================================
% Colormaps
% ============================================================
% if isfield(RT,'rt_colormaps') && isfield(RT.rt_colormaps,'sst2')
    cmap_temp = RT.rt_colormaps.greenbrown;
% else
%     cmap_temp = cmocean('thermal',256);
% end

% if isfield(RT,'rt_colormaps') && isfield(RT.rt_colormaps,'section')
    cmap_salt = RT.rt_colormaps.sst2;
% else
%     cmap_salt = cmocean('haline',256);
% end

panel_letters = arrayfun(@(k) ['(' char('a'+k-1) ')'], 1:nDays, 'UniformOutput', false);

%% ============================================================
% Figure 1: Temperature
%% ============================================================
fig1 = figure('Color','w','Units','pixels','Position',[60 60 1150 1100]);
% tl1 = tiledlayout(fig1, nrow, ncol, 'TileSpacing','compact', 'Padding','compact');
tl1 = tiledlayout(fig1, nrow, ncol, 'TileSpacing','compact', 'Padding','compact');
tl1.Position = [0.055 0.06 0.85 0.90];   % leave room on the right for colorbar

ax_last = [];

for ii = 1:nDays
    ax = nexttile(tl1, ii); hold(ax,'on'); box(ax,'on')
    ax_last = ax;

    if isempty(T_all{ii})
        text(ax, 0.5, 0.5, 'No file', 'Units','normalized', ...
            'HorizontalAlignment','center', 'FontSize',12)
        axis(ax, 'off')
        continue
    end

    pcolor(ax, lon_sub, lat_sub, T_all{ii}); shading(ax,'flat')
    colormap(ax, cmap_temp)
    caxis(ax, [Tmin Tmax])

    quiver(ax, ...
        lon_sub(1:skip_x:end), lat_sub(1:skip_y:end), ...
        U_all{ii}(1:skip_y:end,1:skip_x:end), ...
        V_all{ii}(1:skip_y:end,1:skip_x:end), ...
        quiv_scale, 'k');

    plot(ax, ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

    xlim(ax, [lon_min lon_max])
    ylim(ax, [lat_min lat_max])

    % let the axes fill the tile better
    set(ax, 'FontSize',11, 'LineWidth',1, 'TickDir','out')
    set(ax, 'PlotBoxAspectRatio', [1 1 1])

    title(ax, char(time_str(ii)), 'FontSize',12, 'FontWeight','bold')

    text(ax, 0.02, 0.985, panel_letters{ii}, ...
    'Units','normalized', ...
    'FontSize',13, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left');

    if mod(ii-1,ncol) == 0
        ylabel(ax, 'Lat (^oN)', 'FontSize',12)
    else
        set(ax, 'YTickLabel', [])
    end

    if ii > ncol*(nrow-1)
        xlabel(ax, 'Lon (^oE)', 'FontSize',12)
    else
        set(ax, 'XTickLabel', [])
    end
end

% title(tl1, sprintf('SWOT DA daily temperature at -%.0f m, hour = %d, 21 Apr to 6 May 2023', depth_used, hour_id), ...
%     'FontSize',16, 'FontWeight','bold')

cb1 = colorbar('position',[ 0.91 0.3 0.015 0.4]);
% cb1.Layout.Tile = 'east';
cb1.Label.String = sprintf('SWOT DA Temperature at z = -%.0f m (^oC)', depth_used);
cb1.FontSize = 14;
cb1.TickDirection = 'out';
cb1.YAxisLocation = 'right';
print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/si_Fig1'], '-dpng', '-r300');

%% ============================================================
% Figure 2: Salinity
%% ============================================================
fig2 = figure('Color','w','Units','pixels','Position',[90 90 1150 1100]);
% tl2 = tiledlayout(fig2, nrow, ncol, 'TileSpacing','compact', 'Padding','compact');
tl2 = tiledlayout(fig2, nrow, ncol, 'TileSpacing','compact', 'Padding','compact');
tl2.Position = [0.055 0.06 0.85 0.90];   % leave room on the right for colorbar

ax_last = [];

for ii = 1:nDays
    ax = nexttile(tl2, ii); hold(ax,'on'); box(ax,'on')
    ax_last = ax;

    if isempty(S_all{ii})
        text(ax, 0.5, 0.5, 'No file', 'Units','normalized', ...
            'HorizontalAlignment','center', 'FontSize',12)
        axis(ax, 'off')
        continue
    end

    pcolor(ax, lon_sub, lat_sub, S_all{ii}); shading(ax,'flat')
    colormap(ax, cmap_salt)
    caxis(ax, [Smin Smax])

    quiver(ax, ...
        lon_sub(1:skip_x:end), lat_sub(1:skip_y:end), ...
        U_all{ii}(1:skip_y:end,1:skip_x:end), ...
        V_all{ii}(1:skip_y:end,1:skip_x:end), ...
        quiv_scale, 'k');

    plot(ax, ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

    xlim(ax, [lon_min lon_max])
    ylim(ax, [lat_min lat_max])

    set(ax, 'FontSize',11, 'LineWidth',1, 'TickDir','out')
    set(ax, 'PlotBoxAspectRatio', [1 1 1])

    title(ax, char(time_str(ii)), 'FontSize',12, 'FontWeight','bold')

    text(ax, 0.02, 0.985, panel_letters{ii}, ...
    'Units','normalized', ...
    'FontSize',13, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left');

    if mod(ii-1,ncol) == 0
        ylabel(ax, 'Lat (^oN)', 'FontSize',12)
    else
        set(ax, 'YTickLabel', [])
    end

    if ii > ncol*(nrow-1)
        xlabel(ax, 'Lon (^oE)', 'FontSize',12)
    else
        set(ax, 'XTickLabel', [])
    end
end
% colormap(RT.rt_colormaps.sst)
%
% title(tl2, sprintf('SWOT DA daily salinity at -%.0f m, hour = %d, 21 Apr to 6 May 2023', depth_used, hour_id), ...
%     'FontSize',16, 'FontWeight','bold')

cb2 = colorbar('position',[ 0.91 0.3 0.015 0.4]);
% cb2.Layout.Tile = 'east';
cb2.Label.String = sprintf('SWOT DA Salinity at z = -%.0f m', depth_used);
cb2.FontSize = 14;
cb2.TickDirection = 'out';
cb2.YAxisLocation = 'right';
print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/si_Fig2'], '-dpng', '-r300');
