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


clear all; close all; clc

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
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/', ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'};

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


% ind{1} = find(data_Leg.Time > datenum(2023,05,08,02,00,00) &... % test
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
depth_lim_plot = [-205 -40];
vel_lim        = [-30 30];   % cm/s

% For metrics, only use upper 200 m
iz_adcp_metric = find(double(data_Leg.depth) <= 200);

% If you want 50-200 m instead, use:
% iz_adcp_metric = find(double(data_Leg.depth) >= 50 & double(data_Leg.depth) <= 200);




for itr = 1:2

  if itr == 1
        xlim1 = 39.506;
        xlim2 = 39.95;
    elseif itr == 2
        xlim1 = 39.494;
        xlim2 = 39.922;
  end

    JP_ind = ind{itr};

    %% --------------------------------------------------------
    % Transect geometry
    %% --------------------------------------------------------
    latinit = data_Leg.lat(JP_ind(1));
    loninit = data_Leg.lon(JP_ind(1));
    latend  = data_Leg.lat(JP_ind(end));
    lonend  = data_Leg.lon(JP_ind(end));

    dy_sec_km = sign(latend - latinit) * ...
        geodistance(latend, mean([lonend loninit]), latinit, mean([lonend loninit])) / 1000;

    dx_sec_km = sign(lonend - loninit) * ...
        geodistance(mean([latinit latend]), lonend, mean([latinit latend]), loninit) / 1000;

    angle_sec = atan2d(dy_sec_km, dx_sec_km);

    %% --------------------------------------------------------
    % ADCP rotated cross-section velocity
    %% --------------------------------------------------------
    uv_complex_rotated_adcp = ...
        (data_Leg.u(:,JP_ind) + 1i*data_Leg.v(:,JP_ind)) .* exp(1i*(-angle_sec)*pi/180);

    if itr == 1
        uv_crosssection_adcp = -imag(uv_complex_rotated_adcp);   % same convention as before
    elseif itr == 2
        uv_crosssection_adcp = imag(uv_complex_rotated_adcp);
    end
    %% --------------------------------------------------------
    % Time labels
    %% --------------------------------------------------------
    t0   = data_Leg.Time(JP_ind(1));
    t1   = data_Leg.Time(JP_ind(end));
    tmid = mean([t0 t1]);

    %% --------------------------------------------------------
    % Figure
    %% --------------------------------------------------------
    fig = figure('Color','w','Units','pixels','Position',[100 120 1350 400]);
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

        [~, model] = read_nc_file_struct([model_paths{kk} best_name]);
        model.time = datenum(seconds(model.ocean_time) + datenum(1968,5,23));

        clear u_model_at_ADCP v_model_at_ADCP
        for iz = 1:length(model.depth)
            u_model_at_ADCP(:,iz) = interp2(double(model.lon), double(model.lat), ...
                squeeze(model.u(:,:,iz,best_hr))', ...
                data_Leg.lon(JP_ind), data_Leg.lat(JP_ind), 'linear');

            v_model_at_ADCP(:,iz) = interp2(double(model.lon), double(model.lat), ...
                squeeze(model.v(:,:,iz,best_hr))', ...
                data_Leg.lon(JP_ind), data_Leg.lat(JP_ind), 'linear');
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

        RMSE_now(kk) = sqrt(mean((A(valid) - B(valid)).^2)) * 100;   % cm/s
        CORR_now(kk) = corr(A(valid), B(valid), 'rows', 'complete');


        % plot model section on model depth grid to show deeper structure
        ax = nexttile(tl,kk); hold(ax,'on'); box(ax,'on')
        pcolor(ax, data_Leg.lat(JP_ind), -double(model.depth), uv_crosssection_model'*100);
        shading(ax,'interp')
        % cmocean('balance')
        colormap(RT.rt_colormaps.redblueclass)

        caxis(ax, vel_lim)
        ylim(ax, depth_lim_plot)
        set(ax,'Layer','top','FontSize',13,'TickDir','out')

        x1 = min(data_Leg.lat(JP_ind));
        x2 = max(data_Leg.lat(JP_ind));
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

        text(ax, 0.02, 0.95, panel_labels{kk}, ...
            'Units','normalized','FontSize',14,'FontWeight','bold','Color','k')

        text(ax, 0.04, 0.08, ...
            {sprintf('RMSE = %.2f cm s^{-1}', RMSE_now(kk)), ...
             sprintf('r = %.2f', CORR_now(kk))}, ...
            'Units','normalized','FontSize',14,'FontWeight','bold', ...
            'Color','k','Margin',2)

        xlim([xlim1 xlim2])
    end

    %% ========================================================
    % Panel 3: ADCP
    %% ========================================================
    ax3 = nexttile(tl,3); hold(ax3,'on'); box(ax3,'on')
    pcolor(ax3, data_Leg.lat(JP_ind), -double(data_Leg.depth), uv_crosssection_adcp*100);
    shading(ax3,'flat')
    colormap(RT.rt_colormaps.redblueclass)

    caxis(ax3, vel_lim)
    ylim(ax3, depth_lim_plot)
    xlim([xlim1 xlim2])
    set(ax3,'Layer','top','FontSize',13,'TickDir','out')

    title(ax3, ['ADCP: ' datestr(t0,'dd-mmm HH:MM') ' -- ' datestr(t1,'dd-mmm HH:MM')], ...
        'FontSize',13,'FontWeight','bold')
    xlabel(ax3,'Lat (^oN)','FontSize',13)
    set(ax3,'YTickLabel',[])

    text(ax3, 0.02, 0.95, panel_labels{3}, ...
        'Units','normalized','FontSize',14,'FontWeight','bold','Color','k')

    cb = colorbar(ax3,'eastoutside');
    ylabel(cb,'Cross-section velocity (cm s^{-1})','FontSize',12)
    cb.FontSize = 12;
    cb.TickDirection = 'out';

    % %% --------------------------------------------------------
    % % Console output
    % %% --------------------------------------------------------
    % disp(' ')
    % disp(['Transect ' num2str(itr)])
    % disp(['ADCP window: ' datestr(t0,'dd-mmm-yyyy HH:MM') ' -- ' datestr(t1,'dd-mmm-yyyy HH:MM')])
    % fprintf('  SWOT DA  RMSE = %.2f cm s^-1, r = %.2f\n', RMSE_now(1), CORR_now(1))
    % fprintf('  Joint DA RMSE = %.2f cm s^-1, r = %.2f\n', RMSE_now(2), CORR_now(2))
    if itr == 1
        print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_ADCP_SI', '-dpng', '-r300');
    end
    if itr == 2
        print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_ADCP', '-dpng', '-r300');
    end
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