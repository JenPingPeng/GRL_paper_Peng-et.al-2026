% ============================================================
% Figure 4 / Supporting Information release script
% Purpose: Validate forecast-phase simulations against two independent glider sections.
% Outputs: Profile-matched RMSE statistics and section comparison figures.
%
% Notes for release:
% - Run this script directly after confirming that the local data paths are available.
% - Local absolute paths are intentionally preserved from the submitted analysis version.
% - Keep figure settings unchanged unless you are intentionally updating the submitted figure.
% ============================================================

clearvars; close all; clc

%% ============================================================
% Toolboxes and colormaps
%% ============================================================
addpath TEOS_10
addpath TEOS_10/library/
addpath Sea-Bird-Toolbox-master/CTD_CNV/
addpath cmocean

RT = load("rt_colormaps.mat");

set(groot,'defaultAxesFontName','Helvetica')
set(groot,'defaultTextFontName','Helvetica')

%% ============================================================
% Paths
%% ============================================================
glider_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/';

WMOP_paths = { ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v6_fr_07_05/', ... % GEN
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/', ... % SWOT DA
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8_fr_07_05/', ... % In-situ DA
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'  ... % Joint DA
    };

version_labels = {'GEN','SWOT DA','In-situ DA','Joint DA'};

dir_WMOPs = cell(1,4);
for kk = 1:4
    dir_WMOPs{kk} = dir(fullfile(WMOP_paths{kk}, '*.nc'));
end

%% ============================================================
% Build model time catalogs once
%% ============================================================
fprintf('\nBuilding model time catalogs ...\n')
catalogs = cell(1,4);
for kk = 1:4
    catalogs{kk} = build_model_time_catalog(WMOP_paths{kk}, dir_WMOPs{kk});
    fprintf('  %s: %d snapshots indexed\n', version_labels{kk}, numel(catalogs{kk}.time))
end

%% ============================================================
% Load glider data
%% ============================================================
[~, glider1] = read_nc_file_struct([glider_path 'dep0040_sdeep01_scb-sldeep001_L2_2023-04-25_data_dt.nc']);
glider1.Time = datenum(seconds(glider1.time) + datenum(1970,01,01));
glider1.CT = gsw_CT_from_t(glider1.salinity, glider1.temperature, 0);
glider1.SA = gsw_SA_from_SP(glider1.salinity, 0, glider1.longitude, glider1.latitude);
glider1.PD = gsw_rho(glider1.SA, glider1.CT, 0) - 1000;

[~, glider2] = read_nc_file_struct([glider_path 'dep0007_sdeep09_scb-sldeep009_L2_2023-04-25_data_dt.nc']);
glider2.Time = datenum(seconds(glider2.time) + datenum(1970,01,01));
glider2.CT = gsw_CT_from_t(glider2.salinity, glider2.temperature, 0);
glider2.SA = gsw_SA_from_SP(glider2.salinity, 0, glider2.longitude, glider2.latitude);
glider2.PD = gsw_rho(glider2.SA, glider2.CT, 0) - 1000;

%% ============================================================
% Define glider sections
%% ============================================================
g1_ind = find(glider1.Time >= datenum(2023,5,7,23,40,00) & ...
              glider1.Time <= datenum(2023,5,10,07,03,00));

g2_ind = find(glider2.Time >= datenum(2023,5,7,07,38,00) & ...
              glider2.Time <= datenum(2023,5,10,07,12,00));

% g1_ind = find(glider1.Time >= datenum(2023,5,7,23,40,00) & ...
%               glider1.Time <= datenum(2023,5,10,0,00,00));
%
% g2_ind = find(glider2.Time >= datenum(2023,5,7,07,38,00) & ...
%               glider2.Time <= datenum(2023,5,10,00,00,00));
%% ============================================================
% Plot / contour settings
%% ============================================================
lat_xlim   = [39.54 40.15];
depth_ylim = [-700 -0];

temp_lim   = [13.1 13.8];
salt_lim   = [38.2 38.6];
rho_lim    = [28.80 29.10];

rho_cont_all   = 28.8:0.05:29.1;
rho_cont_label = [28.9 28.95];

if isfield(RT,'rt_colormaps') && isfield(RT.rt_colormaps,'sst2')
    cmap_temp = RT.rt_colormaps.sst2;
else
    cmap_temp = cmocean('thermal',256);
end

if isfield(RT,'rt_colormaps') && isfield(RT.rt_colormaps,'section')
    cmap_salt = RT.rt_colormaps.section;
else
    cmap_salt = parula(256);
end

cmap_rho = cmocean('dense',256);

%% ============================================================
% Process each glider using profile-matched model times
%% ============================================================
out_g1 = process_glider_section_profilematched( ...
    glider1, g1_ind, 'Glider 1', ...
    WMOP_paths, dir_WMOPs, catalogs, version_labels, ...
    lat_xlim, depth_ylim, temp_lim, salt_lim, rho_lim, ...
    cmap_temp, cmap_salt, cmap_rho, ...
    rho_cont_all, rho_cont_label);

out_g2 = process_glider_section_profilematched( ...
    glider2, g2_ind, 'Glider 2', ...
    WMOP_paths, dir_WMOPs, catalogs, version_labels, ...
    lat_xlim, depth_ylim, temp_lim, salt_lim, rho_lim, ...
    cmap_temp, cmap_salt, cmap_rho, ...
    rho_cont_all, rho_cont_label);

%% ============================================================
% Summary table: mean RMSE ± std across profiles
%% ============================================================
all_prof_T  = [out_g1.profRMSE_T;  out_g2.profRMSE_T];
all_prof_S  = [out_g1.profRMSE_S;  out_g2.profRMSE_S];
all_prof_PD = [out_g1.profRMSE_PD; out_g2.profRMSE_PD];

% summary_table = table( ...
%     version_labels(:), ...
%     out_g1.meanRMSE_T(:),  out_g1.stdRMSE_T(:),  ...
%     out_g2.meanRMSE_T(:),  out_g2.stdRMSE_T(:),  ...
%     nanmean(all_prof_T,1)',  nanstd(all_prof_T,0,1)', ...
%     out_g1.meanRMSE_S(:),  out_g1.stdRMSE_S(:),  ...
%     out_g2.meanRMSE_S(:),  out_g2.stdRMSE_S(:),  ...
%     nanmean(all_prof_S,1)',  nanstd(all_prof_S,0,1)', ...
%     out_g1.meanRMSE_PD(:), out_g1.stdRMSE_PD(:), ...
%     out_g2.meanRMSE_PD(:), out_g2.stdRMSE_PD(:), ...
%     nanmean(all_prof_PD,1)', nanstd(all_prof_PD,0,1)', ...
%     'VariableNames', { ...
%     'Simulation', ...
%     'Temp_Glider1_Mean', 'Temp_Glider1_Std', ...
%     'Temp_Glider2_Mean', 'Temp_Glider2_Std', ...
%     'Temp_All_Mean',    'Temp_All_Std', ...
%     'Salt_Glider1_Mean', 'Salt_Glider1_Std', ...
%     'Salt_Glider2_Mean', 'Salt_Glider2_Std', ...
%     'Salt_All_Mean',    'Salt_All_Std', ...
%     'Density_Glider1_Mean', 'Density_Glider1_Std', ...
%     'Density_Glider2_Mean', 'Density_Glider2_Std', ...
%     'Density_All_Mean',    'Density_All_Std'});

summary_table = table( ...
    version_labels(:), ...
    nanmean(all_prof_T,1)',  nanstd(all_prof_T,0,1)', ...
    nanmean(all_prof_S,1)',  nanstd(all_prof_S,0,1)', ...
    nanmean(all_prof_PD,1)', nanstd(all_prof_PD,0,1)', ...
    'VariableNames', { ...
    'Simulation', ...
    'Temp_All_Mean',    'Temp_All_Std', ...
    'Salt_All_Mean',    'Salt_All_Std', ...
    'Density_All_Mean',    'Density_All_Std'});

disp(' ')
disp('Profile-matched RMSE summary across two gliders:')
disp(summary_table)

% optional export
% writetable(summary_table,'Glider_profile_matched_RMSE_summary.csv')

%% ============================================================
% Optional: profile time-offset summary
%% ============================================================
fprintf('\nMean absolute time offset between model and glider profiles:\n')
for kk = 1:4
    fprintf('  %-10s  G1 = %.2f h,  G2 = %.2f h\n', version_labels{kk}, ...
        nanmean(abs(out_g1.dt_hours(:,kk))), nanmean(abs(out_g2.dt_hours(:,kk))))
end

%% ============================================================
% Local functions
%% ============================================================

function out = process_glider_section_profilematched( ...
    glider, ind_sec, glider_name, ...
    WMOP_paths, dir_WMOPs, catalogs, version_labels, ...
    lat_xlim, depth_ylim, temp_lim, salt_lim, rho_lim, ...
    cmap_temp, cmap_salt, cmap_rho, ...
    rho_cont_all, rho_cont_label)

    ind_sec = ind_sec(:);
    nSim = numel(version_labels);

    %% --------------------------------------------------------
    % Section info
    %% --------------------------------------------------------
    lat_tran  = double(glider.latitude(ind_sec));
    lon_tran  = double(glider.longitude(ind_sec));
    time_tran = double(glider.Time(ind_sec));

    glider_t0 = time_tran(1);
    glider_t1 = time_tran(end);

    %% --------------------------------------------------------
    % Get model depth from first model file
    %% --------------------------------------------------------
    test_file = fullfile(WMOP_paths{1}, dir_WMOPs{1}(1).name);
    model_depth = double(ncread(test_file,'depth'));
    nZ = numel(model_depth);
    nP = numel(ind_sec);

    %% --------------------------------------------------------
    % Allocate
    %% --------------------------------------------------------
    T_model  = cell(1,nSim);
    S_model  = cell(1,nSim);
    PD_model = cell(1,nSim);

    matched_times = nan(nP,nSim);
    dt_hours      = nan(nP,nSim);

    for kk = 1:nSim
        T_model{kk}  = nan(nP,nZ);
        S_model{kk}  = nan(nP,nZ);
        PD_model{kk} = nan(nP,nZ);
    end

    %% --------------------------------------------------------
    % Interpolate glider onto model depth for each profile
    %% --------------------------------------------------------
    T_glider_on_model  = nan(nP,nZ);
    S_glider_on_model  = nan(nP,nZ);
    PD_glider_on_model = nan(nP,nZ);

    for ip = 1:nP
        T_glider_on_model(ip,:) = interp1( ...
            double(glider.depth(:)), ...
            double(glider.temperature(:,ind_sec(ip))), ...
            model_depth, 'linear', NaN);

        S_glider_on_model(ip,:) = interp1( ...
            double(glider.depth(:)), ...
            double(glider.salinity(:,ind_sec(ip))), ...
            model_depth, 'linear', NaN);

        PD_glider_on_model(ip,:) = interp1( ...
            double(glider.depth(:)), ...
            double(glider.PD(:,ind_sec(ip))), ...
            model_depth, 'linear', NaN);
    end

    %% --------------------------------------------------------
    % Sample model at nearest time to each glider profile
    %% --------------------------------------------------------
    for kk = 1:nSim
        fprintf('\n%s | %s: sampling profile-matched model fields ...\n', glider_name, version_labels{kk})

        last_file_idx = NaN;
        model_now = [];

        for ip = 1:nP
            target_time = time_tran(ip);

            [file_idx, time_idx, best_time] = find_nearest_snapshot_in_catalog(catalogs{kk}, target_time);

            matched_times(ip,kk) = best_time;
            dt_hours(ip,kk) = 24 * (best_time - target_time);

            if isnan(last_file_idx) || file_idx ~= last_file_idx
                fname = fullfile(WMOP_paths{kk}, dir_WMOPs{kk}(file_idx).name);
                [~, model_now] = read_nc_file_struct(fname);
                model_now.time = datenum(seconds(model_now.ocean_time) + datenum(1968,5,23));
                model_now.rho  = gsw_rho(model_now.salt, model_now.temp, 40) - 1000;
                last_file_idx = file_idx;
            end

            lon0 = lon_tran(ip);
            lat0 = lat_tran(ip);

            for iz = 1:nZ
                T_model{kk}(ip,iz) = interp2( ...
                    double(model_now.lon), double(model_now.lat), ...
                    squeeze(model_now.temp(:,:,iz,time_idx))', ...
                    lon0, lat0, 'linear');

                S_model{kk}(ip,iz) = interp2( ...
                    double(model_now.lon), double(model_now.lat), ...
                    squeeze(model_now.salt(:,:,iz,time_idx))', ...
                    lon0, lat0, 'linear');

                PD_model{kk}(ip,iz) = interp2( ...
                    double(model_now.lon), double(model_now.lat), ...
                    squeeze(model_now.rho(:,:,iz,time_idx))', ...
                    lon0, lat0, 'linear');
            end
        end
    end

    %% --------------------------------------------------------
    % Profile-wise RMSE over 0-600 m
    %% --------------------------------------------------------
    iz_eval = find(model_depth >= 0 & model_depth <= 600);

    profRMSE_T  = nan(nP,nSim);
    profRMSE_S  = nan(nP,nSim);
    profRMSE_PD = nan(nP,nSim);

    for kk = 1:nSim
        for ip = 1:nP
            profRMSE_T(ip,kk) = calc_rmse( ...
                T_model{kk}(ip,iz_eval), T_glider_on_model(ip,iz_eval));

            profRMSE_S(ip,kk) = calc_rmse( ...
                S_model{kk}(ip,iz_eval), S_glider_on_model(ip,iz_eval));

            profRMSE_PD(ip,kk) = calc_rmse( ...
                PD_model{kk}(ip,iz_eval), PD_glider_on_model(ip,iz_eval));
        end
    end

    meanRMSE_T  = nanmean(profRMSE_T,1);
    stdRMSE_T   = nanstd(profRMSE_T,0,1);

    meanRMSE_S  = nanmean(profRMSE_S,1);
    stdRMSE_S   = nanstd(profRMSE_S,0,1);

    meanRMSE_PD = nanmean(profRMSE_PD,1);
    stdRMSE_PD  = nanstd(profRMSE_PD,0,1);

    %% --------------------------------------------------------
    % Real glider fields on native depth grid
    %% --------------------------------------------------------
    [lat_unique, idx_unique] = unique(lat_tran, 'stable');
    lat_glider_real = lat_unique;
    ind_unique = ind_sec(idx_unique);

    T_glider_real  = double(glider.temperature(:,ind_unique));
    S_glider_real  = double(glider.salinity(:,ind_unique));
    PD_glider_real = double(glider.PD(:,ind_unique));

    %% --------------------------------------------------------
    % Make figures
    %% --------------------------------------------------------
    make_section_figure( ...
        T_model, T_glider_real, PD_model, PD_glider_real, ...
        lat_tran, lat_glider_real, model_depth, double(glider.depth(:)), ...
        matched_times, dt_hours, ...
        version_labels, glider_name, glider_t0, glider_t1, ...
        meanRMSE_T, stdRMSE_T, ...
        'Temperature', '°C', temp_lim, cmap_temp, ...
        lat_xlim, depth_ylim, rho_cont_all, rho_cont_label, true);

    make_section_figure( ...
        S_model, S_glider_real, PD_model, PD_glider_real, ...
        lat_tran, lat_glider_real, model_depth, double(glider.depth(:)), ...
        matched_times, dt_hours, ...
        version_labels, glider_name, glider_t0, glider_t1, ...
        meanRMSE_S, stdRMSE_S, ...
        'Salinity', 'g kg^{-1}', salt_lim, cmap_salt, ...
        lat_xlim, depth_ylim, rho_cont_all, rho_cont_label, true);

    make_section_figure( ...
        PD_model, PD_glider_real, PD_model, PD_glider_real, ...
        lat_tran, lat_glider_real, model_depth, double(glider.depth(:)), ...
        matched_times, dt_hours, ...
        version_labels, glider_name, glider_t0, glider_t1, ...
        meanRMSE_PD, stdRMSE_PD, ...
        'Density', 'kg m^{-3}', rho_lim, cmap_rho, ...
        lat_xlim, depth_ylim, rho_cont_all, rho_cont_label, false);

    %% --------------------------------------------------------
    % Return outputs
    %% --------------------------------------------------------
    out.profRMSE_T  = profRMSE_T;
    out.profRMSE_S  = profRMSE_S;
    out.profRMSE_PD = profRMSE_PD;

    out.meanRMSE_T  = meanRMSE_T;
    out.stdRMSE_T   = stdRMSE_T;

    out.meanRMSE_S  = meanRMSE_S;
    out.stdRMSE_S   = stdRMSE_S;

    out.meanRMSE_PD = meanRMSE_PD;
    out.stdRMSE_PD  = stdRMSE_PD;

    out.dt_hours = dt_hours;
    out.matched_times = matched_times;
end

function make_section_figure( ...
    model_cell, glider_real, PD_model_cell, PD_glider_real, ...
    lat_tran, lat_glider_real, model_depth, glider_depth, ...
    matched_times, dt_hours, ...
    version_labels, glider_name, glider_t0, glider_t1, ...
    meanRMSE, stdRMSE, ...
    var_name, var_unit, clim_now, cmap_now, ...
    lat_xlim, depth_ylim, rho_cont_all, rho_cont_label, do_contours)

    nSim = numel(version_labels);

    fig = figure('Color','w','Units','pixels','Position',[80 120 2800 520]);
    tl  = tiledlayout(fig,1,nSim+1,'TileSpacing','compact','Padding','compact');

    panel_letters = {'(a)','(b)','(c)','(d)','(e)'};

    for kk = 1:(nSim+1)
        ax = nexttile(tl,kk);
        hold(ax,'on')
        box(ax,'on')

        if kk <= nSim
            pcolor(ax, lat_tran, -model_depth, model_cell{kk}');
            shading(ax,'flat')
            colormap(ax, cmap_now)
            caxis(ax, clim_now)

            if do_contours
                contour(ax, lat_tran, -model_depth, PD_model_cell{kk}', ...
                    rho_cont_all, 'k', 'LineWidth', 0.45);
                [C,h] = contour(ax, lat_tran, -model_depth, PD_model_cell{kk}', ...
                    rho_cont_label, 'k', 'LineWidth', 0.9);
                clabel(C,h,'FontSize',11,'Color','k','LabelSpacing',350)
            end

            title(ax, sprintf('%s\nmean|dt| = %.2f h', ...
                version_labels{kk}, nanmean(abs(dt_hours(:,kk)))), ...
                'FontSize',14,'FontWeight','bold')

            text(ax, 0.02, 0.95, panel_letters{kk}, 'Units','normalized', ...
                'FontSize',14,'FontWeight','bold','Color','k')

            text(ax, 0.04, 0.06, sprintf('RMSE = %.3f ± %.3f %s', ...
                meanRMSE(kk), stdRMSE(kk), var_unit), ...
                'Units','normalized', 'FontSize',13, 'FontWeight','bold', ...
                'Color','k', 'BackgroundColor','w', 'Margin',2)

        else
            pcolor(ax, lat_glider_real, -glider_depth, glider_real);
            shading(ax,'flat')
            colormap(ax, cmap_now)
            caxis(ax, clim_now)

            if do_contours
                contour(ax, lat_glider_real, -glider_depth, PD_glider_real, ...
                    rho_cont_all, 'k', 'LineWidth', 0.45);
                [C,h] = contour(ax, lat_glider_real, -glider_depth, PD_glider_real, ...
                    rho_cont_label, 'k', 'LineWidth', 0.9);
                clabel(C,h,'FontSize',11,'Color','k','LabelSpacing',350)
            end

            title(ax, sprintf('%s\n%s to %s', ...
                glider_name, datestr(glider_t0,'dd-mmm HH:MM'), datestr(glider_t1,'dd-mmm HH:MM')), ...
                'FontSize',14,'FontWeight','bold')

            text(ax, 0.02, 0.95, panel_letters{kk}, 'Units','normalized', ...
                'FontSize',14,'FontWeight','bold','Color','k')
        end

        xlim(ax, lat_xlim)
        ylim(ax, depth_ylim)
        set(ax,'FontSize',13,'LineWidth',1,'TickDir','out','Layer','top')

        xlabel(ax,'Lat (\circN)','FontSize',14)
        if kk == 1
            ylabel(ax,'Depth (m)','FontSize',14)
        else
            set(ax,'YTickLabel',[])
        end
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = sprintf('%s (%s)', var_name, var_unit);
    cb.FontSize = 12;
    cb.TickDirection = 'out';

    sgtitle(fig, sprintf('%s | profile-matched comparison', var_name), ...
        'FontSize',16, 'FontWeight','bold')
end

function rmse_val = calc_rmse(A, B)
    mask = isfinite(A) & isfinite(B);
    if nnz(mask) < 3
        rmse_val = NaN;
    else
        d = A(mask) - B(mask);
        rmse_val = sqrt(mean(d.^2));
    end
end

function catalog = build_model_time_catalog(model_path, dir_model)
    all_time = [];
    all_file_idx = [];
    all_time_idx = [];

    for ii = 1:length(dir_model)
        fname = fullfile(model_path, dir_model(ii).name);

        try
            ocean_time = ncread(fname,'ocean_time');
        catch
            continue
        end

        tt = datenum(seconds(ocean_time) + datenum(1968,5,23));

        ntt = numel(tt);
        all_time     = [all_time; tt(:)];
        all_file_idx = [all_file_idx; repmat(ii, ntt, 1)];
        all_time_idx = [all_time_idx; (1:ntt)'];
    end

    [all_time, isort] = sort(all_time);
    all_file_idx = all_file_idx(isort);
    all_time_idx = all_time_idx(isort);

    catalog.time = all_time;
    catalog.file_idx = all_file_idx;
    catalog.time_idx = all_time_idx;
end

function [file_idx, time_idx, best_time] = find_nearest_snapshot_in_catalog(catalog, target_time)
    [~, i0] = min(abs(catalog.time - target_time));
    file_idx = catalog.file_idx(i0);
    time_idx = catalog.time_idx(i0);
    best_time = catalog.time(i0);
end
