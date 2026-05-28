% ============================================================
% Figure 4 release script
% Panels: Figure 4d-f
% Purpose: Cross-eddy velocity transects from SWOT DA, Joint DA, and ADCP.
%
% Notes for release:
% - Run this script for the ADCP validation panels.
% - This script corresponds only to panels (d-f).
% - It computes RMSE and correlation against ADCP for the selected transects.
% ============================================================

clearvars; close all; clc

%% ============================================================
% ADCP validation for two later transects
% Compare: ADCP, SWOT DA, Joint DA
%% ============================================================

addpath TEOS_10
addpath TEOS_10/library/
addpath Sea-Bird-Toolbox-master/CTD_CNV/
addpath cmocean

RT = load("rt_colormaps.mat");

colors =[
          0.2310    0.6660    0.1960

          0.9608    0.4667    0.1608 ];

%% ============================================================
% Paths
%% ============================================================
adcp_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/ADCP/';

model_names = {'SWOT DA','Joint DA'};
model_paths = { ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8_fr_07_05/', ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'};
% model_paths = { ...
%     '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/', ...
%     '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v11_fr_07_05/'};

% panel_labels = {'(a)','(b)','(c)'};

panel_labels = {'(d)','(e)','(f)','(d)','(e)','(f)'};

%% ============================================================
% Load ADCP
%% ============================================================
[~, data_Leg] = read_nc_file_struct([adcp_path 'ADCP_OS150_FaSt-SWOT_LEG2.nc']);
data_Leg.Time = data_Leg.time + datenum(2023,1,1);

%% ============================================================
% Define the two target transects
%% ============================================================
ind{1} = find(data_Leg.Time > datenum(2023,05,08,20,30,00) & ...
              data_Leg.Time < datenum(2023,05,09,03,35,00));

ind{2} = find(data_Leg.Time > datenum(2023,05,09,04,40,00) & ...
              data_Leg.Time < datenum(2023,05,09,14,30,00));

% ind{1} = find(data_Leg.Time > datenum(2023,05,08,02,00,00) &...
%               data_Leg.Time < datenum(2023,05,08,09,35,00));
% %
% ind{2} = find(data_Leg.Time > datenum(2023,05,08,11,00,00) &...
%               data_Leg.Time < datenum(2023,05,08,19,35,00));
%
% ind{3} = find(data_Leg.Time > datenum(2023,05,08,20,30,00) &...
%               data_Leg.Time < datenum(2023,05,09,03,35,00));
%
% ind{4} = find(data_Leg.Time > datenum(2023,05,09,04,40,00) &...
%               data_Leg.Time < datenum(2023,05,09,14,30,00));

% ind{2} = find(data_Leg.Time > datenum(2023,05,09,15,40,00) &...
%               data_Leg.Time < datenum(2023,05,09,23,30,00));

%% ============================================================
% Plot settings
%% ============================================================
depth_lim_plot = [-700 -30];
vel_lim        = [-30 30];   % cm/s

% For metrics, only use upper 200 m
iz_adcp_metric = find(double(data_Leg.depth) <= 200);

% If you want 50-200 m instead, use:
% iz_adcp_metric = find(double(data_Leg.depth) >= 50 & double(data_Leg.depth) <= 200);

for itr = 1 :2

  if itr == 1
        xlim1 = 39.506;
        xlim2 = 39.95;
    elseif itr == 2
        xlim1 = 39.494;
        xlim2 = 39.922;
  end

    transect_ind = ind{itr};

    %% --------------------------------------------------------
    % Transect geometry
    %% --------------------------------------------------------
    latinit = data_Leg.lat(transect_ind(1));
    loninit = data_Leg.lon(transect_ind(1));
    latend  = data_Leg.lat(transect_ind(end));
    lonend  = data_Leg.lon(transect_ind(end));

    dy_sec_km = sign(latend - latinit) * ...
        geodistance(latend, mean([lonend loninit]), latinit, mean([lonend loninit])) / 1000;

    dx_sec_km = sign(lonend - loninit) * ...
        geodistance(mean([latinit latend]), lonend, mean([latinit latend]), loninit) / 1000;

    angle_sec = atan2d(dy_sec_km, dx_sec_km);

    %% --------------------------------------------------------
    % ADCP rotated cross-section velocity
    %% --------------------------------------------------------
    uv_complex_rotated_adcp = ...
        (data_Leg.u(:,transect_ind) + 1i*data_Leg.v(:,transect_ind)) .* exp(1i*(-angle_sec)*pi/180);

    if itr == 1
        uv_crosssection_adcp = -imag(uv_complex_rotated_adcp);   % same convention as before
    elseif itr == 2
        uv_crosssection_adcp = imag(uv_complex_rotated_adcp);
    end
    %% --------------------------------------------------------
    % Time labels
    %% --------------------------------------------------------
    t0   = data_Leg.Time(transect_ind(1));
    t1   = data_Leg.Time(transect_ind(end));
    tmid = mean([t0 t1]);

    %% --------------------------------------------------------
    % Figure
    %% --------------------------------------------------------
    % fig = figure('Color','w','Units','pixels','Position',[100 120 1400 420]);
    fig = figure('Color','w','Units','pixels','Position',[90 120 1150 420]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    %% ========================================================
    % Panels 1-2: models
    %% ========================================================
    RMSE_now = nan(1,2);
    CORR_now = nan(1,2);

    for kk = 1:2

        % -----------------------------------------------------
        % Find nearest model snapshot automatically
        % -----------------------------------------------------
        [best_file, best_hr, best_time, best_name] = find_nearest_snapshot(model_paths{kk}, tmid);
        if itr == 1
            best_name = 'roms_regular_his_20230509.nc';
            best_hr = 4 ; % 21 is good
            best_time= datenum(2023,5,9,best_hr-1,0,0);
        elseif itr == 2
            best_name = 'roms_regular_his_20230509.nc';
            best_hr = 6; % 6 is good
            best_time= datenum(2023,5,9,best_hr-1,0,0);
        end
        [~, model] = read_nc_file_struct([model_paths{kk} best_name]);
        model.time = datenum(seconds(model.ocean_time) + datenum(1968,5,23));

        clear u_model_at_ADCP v_model_at_ADCP
        for iz = 1:length(model.depth)
            u_model_at_ADCP(:,iz) = interp2(double(model.lon), double(model.lat), ...
                squeeze(model.u(:,:,iz,best_hr))', ...
                data_Leg.lon(transect_ind), data_Leg.lat(transect_ind), 'linear');

            v_model_at_ADCP(:,iz) = interp2(double(model.lon), double(model.lat), ...
                squeeze(model.v(:,:,iz,best_hr))', ...
                data_Leg.lon(transect_ind), data_Leg.lat(transect_ind), 'linear');
        end

        uv_complex_rotated_model = ...
            (u_model_at_ADCP + 1i*v_model_at_ADCP) .* exp(1i*(-angle_sec)*pi/180);

        % uv_crosssection_model = imag(uv_complex_rotated_model);   % nTrack x nDepth

        if itr == 1
            uv_crosssection_model = -imag(uv_complex_rotated_model);   % same convention as before
        elseif itr == 2
            uv_crosssection_model = imag(uv_complex_rotated_model);
        end

        % interpolate model to ADCP depth for metric only
        clear uv_crosssection_model_ADCP
        for ip = 1:size(uv_crosssection_model,1)
            uv_crosssection_model_ADCP(:,ip) = interp1(double(model.depth), ...
                uv_crosssection_model(ip,:), double(data_Leg.depth), 'linear', NaN);
        end

        % metrics only within -205 m to -20 m
        iz_adcp_metric = find(double(data_Leg.depth) >= 40 & double(data_Leg.depth) <= 205);

        A = double(uv_crosssection_model_ADCP(iz_adcp_metric,:));
        B = double(uv_crosssection_adcp(iz_adcp_metric,:));

        valid = isfinite(A) & isfinite(B);

        % RMSE_now(kk) = sqrt(mean((A(valid) - B(valid)).^2)) * 100;   % cm/s
        err = (A(valid) - B(valid)) * 100;   % cm/s

        RMSE_now(kk) = sqrt(mean(err.^2, 'omitnan'));
        STD_now(kk)  = std(err, 'omitnan');
        CORR_now(kk) = corr(A(valid), B(valid), 'rows', 'complete');

        % plot model section on model depth grid to show deeper structure
        ax = nexttile(tl,kk); hold(ax,'on'); box(ax,'on')
        pcolor(ax, data_Leg.lat(transect_ind), -double(model.depth), uv_crosssection_model'*100);
        x = double(data_Leg.lat(transect_ind));
        y = -double(model.depth(:));
        Z = uv_crosssection_model' * 100;   % depth x track

        % Sort by latitude
        [xs, isort] = sort(x);
        Zs = Z(:, isort);

        % Remove repeated latitude values
        [xu, ia] = unique(xs, 'stable');
        Zu = Zs(:, ia);

        [C,h] = contour(ax, xu, y, Zu, -20:5:20, ...
            'k-', 'LineWidth', 0.8);

        clabel(C,h,'FontSize',12,'Color','k')

        shading(ax,'interp')
        % cmocean('balance')
        colormap(RT.rt_colormaps.redblueclass)

        caxis(ax, vel_lim)
        ylim(ax, depth_lim_plot)
        set(ax,'Layer','top','FontSize',13,'TickDir','out','ytick',[-700:100:0])

        x1 = min(data_Leg.lat(transect_ind));
        x2 = max(data_Leg.lat(transect_ind));
        metric_z_top = -30;
        metric_z_bot = -205;
        % rectangle(ax, 'Position', [x1, metric_z_bot, x2-x1, (metric_z_top-metric_z_bot)*0.98], ...
        %     'EdgeColor', colors(kk,:), ...
        %     'LineStyle', '-', ...
        %     'LineWidth', 2);

        title(ax, sprintf('%s: %s', model_names{kk}, datestr(best_time,'dd-mmm HH:MM')), ...
            'FontSize',13,'FontWeight','bold')
        xlabel(ax,'Lat (^oN)','FontSize',13)

        if kk == 1
            ylabel(ax,'Depth (m)','FontSize',13)
        else
            set(ax,'YTickLabel',[])
        end

        text(ax, 0.04, 0.95, panel_labels{kk}, ...
            'Units','normalized','FontSize',14,'FontWeight','bold','Color','k')
% Semi-transparent white patch behind text
xp = [39.51 39.81 39.81 39.51];
yp = [-495-180 -495-180 -420-180 -420-180];

patch(ax, xp, yp, 'w', ...
    'EdgeColor','none', ...
    'FaceAlpha',0.85);

text(ax, 0.04, 0.09, ...
    {sprintf('RMSE = %.2f cm s^{-1}', RMSE_now(kk)), ...
     sprintf('r = %.2f', CORR_now(kk))}, ...
    'Units','normalized', ...
    'FontSize',14, ...
    'FontWeight','bold', ...
    'Color','k');
        % text(ax, 0.04, 0.12, ...
        %     {sprintf('RMSE = %.2f cm s^{-1}', RMSE_now(kk)), sprintf('r = %.2f', CORR_now(kk))},...
        %     'Units','normalized','FontSize',14,'FontWeight','bold', ...
        %     'Color','k','Margin',2,'BackgroundColor','w')
        % %
        % text(ax, 0.51, 0.05, ...
        %     {sprintf('r = %.2f', CORR_now(kk))}, ...
        %     'Units','normalized','FontSize',14,'FontWeight','bold', ...
        %     'Color','k','Margin',2)
        xlim([xlim1 xlim2])
    end

    %% ========================================================
    % Panel 3: ADCP
    %% ========================================================
  % --------------------------------------------------------
% --------------------------------------------------------
% Interpolate ADCP cross-section velocity onto model section grid
% ADCP: uv_crosssection_adcp is depth x track
% Target grid: model latitude grid x model.depth
% --------------------------------------------------------

lat_adcp = double(data_Leg.lat(transect_ind));     % ADCP track latitude
z_adcp   = double(data_Leg.depth(:));        % ADCP depth, positive downward
Uadcp    = double(uv_crosssection_adcp);     % depth x track, m/s

% Sort ADCP section by latitude, required for interp2
[lat_adcp_sort, isort] = sort(lat_adcp);
Uadcp_sort = Uadcp(:,isort);

% Target model horizontal grid along this section
lat_model_sec = double(model.lat);
lat_model_sec = lat_model_sec(lat_model_sec >= min(lat_adcp_sort) & ...
                              lat_model_sec <= max(lat_adcp_sort));

% Target model vertical grid
z_model = double(model.depth(:));

% Build source and target grids
[LAT_adcp, Z_adcp] = meshgrid(lat_adcp_sort, z_adcp);
[LAT_model, Z_model] = meshgrid(lat_model_sec, z_model);

% Interpolate ADCP onto model latitude-depth grid
Uadcp_on_model = interp2(LAT_adcp, Z_adcp, Uadcp_sort, ...
                         LAT_model, Z_model, 'linear', NaN);

% --------------------------------------------------------
% Smooth horizontally only with a 2 km Gaussian window
% --------------------------------------------------------
win_deg = 2 / 111;   % ~2 km in latitude degrees
dlat_model = median(abs(diff(lat_model_sec)), 'omitnan');
win_pts = max(3, round(win_deg / dlat_model));

% Make window length odd
if mod(win_pts,2) == 0
    win_pts = win_pts + 1;
end

Uadcp_on_model_smooth = smoothdata(Uadcp_on_model, 2, ...
    'gaussian', win_pts, 'omitnan');

% --------------------------------------------------------
% Plot ADCP on model depth grid
% --------------------------------------------------------
ax3 = nexttile(tl,3); hold(ax3,'on'); box(ax3,'on')

pcolor(ax3, lat_model_sec, -z_model, Uadcp_on_model_smooth * 100);
[C,h] = contour(ax3, lat_model_sec, -z_model, Uadcp_on_model_smooth * 100, -40:10:40, 'w-', 'LineWidth', 0.8);
clabel(C,h,'FontSize',12,'Color','w')

shading(ax3,'interp')
colormap(ax3, RT.rt_colormaps.redblueclass)

caxis(ax3, vel_lim)
ylim(ax3, depth_lim_plot)
xlim(ax3, [xlim1 xlim2])
set(ax3,'Layer','top','FontSize',13,'TickDir','out','ytick',[-500:100:0])

% figure; hold on
% plot(lat_model_sec, Uadcp_on_model_smooth(18,:),'b')
% plot(data_Leg.lat(transect_ind), uv_crosssection_model(:,18)','k')

% %% original one
%     ax3 = nexttile(tl,3); hold(ax3,'on'); box(ax3,'on')
%     pcolor(ax3, data_Leg.lat(transect_ind), -double(data_Leg.depth), uv_crosssection_adcp*100);
%
%     shading(ax3,'interp')
%     colormap(RT.rt_colormaps.redblueclass)
%
%
%     % --- contour overlay ---
% x = double(data_Leg.lat(transect_ind));
% y = -double(data_Leg.depth(:));
% Z = double(uv_crosssection_adcp) * 100;   % depth x track
%
% % Sort by latitude for contour
% [xs, isort] = sort(x);
% Zs = Z(:, isort);
%
% % Remove repeated latitudes
% [xu, ia] = unique(xs, 'stable');
% Zu = Zs(:, ia);
%
% [C,h] = contour(ax3, xu, y, Zu, -40:10:40, 'w-', 'LineWidth', 0.8);
% clabel(C,h,'FontSize',12,'Color','w')

    %
    % caxis(ax3, vel_lim)
    % ylim(ax3, depth_lim_plot)
    % xlim([xlim1 xlim2])
    % set(ax3,'Layer','top','FontSize',13,'TickDir','out')

    title(ax3, ['ADCP: ' datestr(t0,'dd-mmm HH:MM') ' - ' datestr(t1,'dd-mmm HH:MM')], ...
        'FontSize',13,'FontWeight','bold')
    xlabel(ax3,'Lat (^oN)','FontSize',13)
    set(ax3,'YTickLabel',[])

    text(ax3, 0.02, 0.95, panel_labels{3}, ...
        'Units','normalized','FontSize',14,'FontWeight','bold','Color','k')

    cb = colorbar(ax3);
    cb.Position = [0.705 0.3 0.2 0.03];
    ylabel(cb,'Cross-section velocity (cm s^{-1})','FontSize',12)
    cb.FontSize = 12;
    cb.TickDirection = 'out';
    cb.Orientation = 'horizontal';

    % if itr == 1
    %     print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_ADCP_SI', '-dpng', '-r300');
    % end
    % if itr == 2
    %     print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_ADCP', '-dpng', '-r300');
    % end
end

%% ============================================================
% Local function: find nearest model snapshot
%% ============================================================
function [best_file, best_hr, best_time, best_name] = find_nearest_snapshot(model_path, target_time)

    dir_all = dir(fullfile(model_path, '*.nc'));

    best_dt   = inf;
    best_file = NaN;
    best_hr   = NaN;
    best_time = NaN;
    best_name = '';

    for ii = 1:length(dir_all)
        fname = fullfile(model_path, dir_all(ii).name);

        try
            ocean_time = ncread(fname, 'ocean_time');
        catch
            continue
        end

        tt = datenum(seconds(ocean_time) + datenum(1968,5,23));

        [dt_min, ih] = min(abs(tt - target_time));

        if dt_min < best_dt
            best_dt   = dt_min;
            best_file = ii;
            best_hr   = ih;
            best_time = tt(ih);
            best_name = dir_all(ii).name;
        end
    end
end
