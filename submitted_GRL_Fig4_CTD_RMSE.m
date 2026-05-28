% ============================================================
% Figure 4 / Supporting Information release script
% Purpose: Validate forecast-phase simulations against independent Leg 2 CTD profiles.
% Outputs: Temperature, salinity, and density profile comparisons, plus RMSE tables.
%
% Notes for release:
% - Run this script directly after confirming that the local data paths are available.
% - Local absolute paths are intentionally preserved from the submitted analysis version.
% - Keep figure settings unchanged unless you are intentionally updating the submitted figure.
% ============================================================

clearvars; close all; clc

%% ============================================================
% Model validation against independent measurements - CTD Leg 2
% Compare CTD profiles against v6-v9 _fr_07_05 simulations
% Figures:
%   1. Temperature profiles (3x3)
%   2. Salinity profiles (3x3)
%   3. Density profiles (3x3)
% Tables:
%   1. Per-profile RMSE
%   2. Mean and std across 9 casts
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
CTD_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/CTDs/Leg2/';

WMOP_paths = { ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v6_fr_07_05/', ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/', ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8_fr_07_05/', ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'};

model_labels = {'GEN','SWOT DA','In-situ DA','Joint DA'};
model_colors = [ ...
    0.0660 0.4430 0.7450
    0.2310 0.6660 0.1960
    0.5210 0.0860 0.8190
    0.9608 0.4667 0.1608];

nSim = numel(WMOP_paths);

dir_models = cell(1,nSim);
for kk = 1:nSim
    dir_models{kk} = dir(fullfile(WMOP_paths{kk}, '*.nc'));
end

% %% ============================================================
% % CTD file list
% %% ============================================================
% CTD_files = dir(fullfile(CTD_path, '*.cnv'));
% if isempty(CTD_files)
%     tmp = dir(CTD_path);
%     CTD_files = tmp(~[tmp.isdir]);
% end
%
% [~,isort] = sort({CTD_files.name});
% CTD_files = CTD_files(isort);
%
% nCast = numel(CTD_files);
% fprintf('Found %d CTD files in Leg 2.\n', nCast)

%% ============================================================
% CTD file list
% Keep only 9 Leg 2 CTD profiles
%% ============================================================
CTD_files = dir(fullfile(CTD_path, '*.cnv'));
if isempty(CTD_files)
    tmp = dir(CTD_path);
    CTD_files = tmp(~[tmp.isdir]);
end

[~,isort] = sort({CTD_files.name});
CTD_files = CTD_files(isort);

% keep only the first 9 casts for Leg 2
if numel(CTD_files) < 9
    error('Expected 9 CTD files for Leg 2, but found only %d.', numel(CTD_files))
elseif numel(CTD_files) > 9
    warning('Found %d CTD files. Using only the first 9 for Leg 2.', numel(CTD_files))
    CTD_files = CTD_files(1:9);
end

nCast = numel(CTD_files);
fprintf('Using %d CTD casts in Leg 2.\n', nCast)

% if nCast ~= 9
%     warning('Expected 9 Leg 2 CTD profiles, but found %d files.', nCast)
% end

%% ============================================================
% Settings
%% ============================================================
depth_max_plot = 700;   % plot to 600 m
depth_max_rmse = 700;   % compute RMSE over overlapping data above 600 m

%% ============================================================
% Containers
%% ============================================================
cast_id   = (1:nCast).';
cast_time = nan(nCast,1);
cast_lat  = nan(nCast,1);
cast_lon  = nan(nCast,1);

obs_depth_all = cell(nCast,1);
obs_T_all     = cell(nCast,1);
obs_S_all     = cell(nCast,1);
obs_PD_all    = cell(nCast,1);

mod_depth_all = cell(nCast,nSim);
mod_T_all     = cell(nCast,nSim);
mod_S_all     = cell(nCast,nSim);
mod_PD_all    = cell(nCast,nSim);
mod_time_all  = nan(nCast,nSim);

RMSE_T  = nan(nCast,nSim);
RMSE_S  = nan(nCast,nSim);
RMSE_PD = nan(nCast,nSim);

allT  = [];
allS  = [];
allPD = [];

%% ============================================================
% Main loop over CTD casts
%% ============================================================
for ic = 1:nCast

    cnv = readSBScnv(fullfile(CTD_path, CTD_files(ic).name));

    [lat_cast, lon_cast, t_num, t_str] = parse_ctd_header(cnv);
    depth_obs = get_ctd_depth(cnv);
    temp_obs  = double(cnv.t090C(:));
    salt_obs  = double(cnv.psal0(:));

    p_obs  = gsw_p_from_z(-depth_obs, lat_cast);
    SA_obs = gsw_SA_from_SP(salt_obs, p_obs, lon_cast, lat_cast);
    CT_obs = gsw_CT_from_t(SA_obs, temp_obs, p_obs);
    pd_obs = gsw_rho(SA_obs, CT_obs, 0) - 1000;

    valid_obs = isfinite(depth_obs) & isfinite(temp_obs) & isfinite(salt_obs) & isfinite(pd_obs);
    depth_obs = double(depth_obs(valid_obs));
    temp_obs  = temp_obs(valid_obs);
    salt_obs  = salt_obs(valid_obs);
    pd_obs    = pd_obs(valid_obs);

    % sort by depth increasing
    [depth_obs, idxs] = sort(depth_obs);
    temp_obs = temp_obs(idxs);
    salt_obs = salt_obs(idxs);
    pd_obs   = pd_obs(idxs);

    cast_time(ic) = t_num;
    cast_lat(ic)  = lat_cast;
    cast_lon(ic)  = lon_cast;

    obs_depth_all{ic} = depth_obs;
    obs_T_all{ic}     = temp_obs;
    obs_S_all{ic}     = salt_obs;
    obs_PD_all{ic}    = pd_obs;

    allT  = [allT; temp_obs(depth_obs <= depth_max_plot)];
    allS  = [allS; salt_obs(depth_obs <= depth_max_plot)];
    allPD = [allPD; pd_obs(depth_obs <= depth_max_plot)];

    for kk = 1:nSim
        [best_name, best_it, best_time] = find_nearest_model_snapshot(WMOP_paths{kk}, dir_models{kk}, t_num);

        [~, model] = read_nc_file_struct(fullfile(WMOP_paths{kk}, best_name));
        model.time = datenum(seconds(model.ocean_time) + datenum(1968,5,23));
        model.rho  = gsw_rho(model.salt, model.temp, 40) - 1000;
        % model.rho  = gsw_rho(model.salt, model.temp, 0) - 1000;

        [depth_mod, T_mod, S_mod, PD_mod] = extract_model_profile(model, best_it, lon_cast, lat_cast);

        mod_depth_all{ic,kk} = depth_mod;
        mod_T_all{ic,kk}     = T_mod;
        mod_S_all{ic,kk}     = S_mod;
        mod_PD_all{ic,kk}    = PD_mod;
        mod_time_all(ic,kk)  = best_time;

        allT  = [allT; T_mod(depth_mod <= depth_max_plot)];
        allS  = [allS; S_mod(depth_mod <= depth_max_plot)];
        allPD = [allPD; PD_mod(depth_mod <= depth_max_plot)];

        % interpolate model to CTD depth for RMSE
        T_mod_i  = interp1(depth_mod, T_mod,  depth_obs, 'linear', NaN);
        S_mod_i  = interp1(depth_mod, S_mod,  depth_obs, 'linear', NaN);
        PD_mod_i = interp1(depth_mod, PD_mod, depth_obs, 'linear', NaN);

        validT = isfinite(depth_obs) & depth_obs <= depth_max_rmse & isfinite(temp_obs) & isfinite(T_mod_i);
        validS = isfinite(depth_obs) & depth_obs <= depth_max_rmse & isfinite(salt_obs) & isfinite(S_mod_i);
        validD = isfinite(depth_obs) & depth_obs <= depth_max_rmse & isfinite(pd_obs)   & isfinite(PD_mod_i);

        if any(validT)
            RMSE_T(ic,kk) = sqrt(mean((T_mod_i(validT) - temp_obs(validT)).^2));
        end
        if any(validS)
            RMSE_S(ic,kk) = sqrt(mean((S_mod_i(validS) - salt_obs(validS)).^2));
        end
        if any(validD)
            RMSE_PD(ic,kk) = sqrt(mean((PD_mod_i(validD) - pd_obs(validD)).^2));
        end
    end
end

%% ============================================================
% Common x-limits for profile figures
%% ============================================================
temp_xlim = [min(allT,[],'omitnan')-0.1,   max(allT,[],'omitnan')+0.1];
salt_xlim = [min(allS,[],'omitnan')-0.02,  max(allS,[],'omitnan')+0.02];
pd_xlim   = [min(allPD,[],'omitnan')-0.05, max(allPD,[],'omitnan')+0.05];

%% ============================================================
% Figure 1: Temperature profiles
%% ============================================================
fig1 = figure('Color','w','Units','pixels','Position',[50 60 1500 1050]);
tl1 = tiledlayout(fig1,3,3,'TileSpacing','compact','Padding','compact');

for ic = 1:nCast
    ax = nexttile(tl1,ic); hold(ax,'on'); box(ax,'on')

    plot(ax, obs_T_all{ic}, -obs_depth_all{ic}, 'k-', 'LineWidth', 2.3)

    for kk = 1:nSim
        plot(ax, mod_T_all{ic,kk}, -mod_depth_all{ic,kk}, '-', ...
            'Color', model_colors(kk,:), 'LineWidth', 1.5)
    end

    xlim(ax, temp_xlim)
    ylim(ax, [-depth_max_plot 0])
    set(ax,'FontSize',11,'LineWidth',1,'TickDir','out')

    title(ax, sprintf('Cast %d  %s', ic, datestr(cast_time(ic),'dd-mmm HH:MM')), ...
        'FontSize',12,'FontWeight','bold')

    if mod(ic-1,3) == 0
        ylabel(ax,'Depth (m)','FontSize',12)
    else
        set(ax,'YTickLabel',[])
    end

    if ic > 6
        xlabel(ax,'Temperature (^oC)','FontSize',12)
    else
        set(ax,'XTickLabel',[])
    end

    if ic == 1
        legend(ax, [{'Obs'}, model_labels], 'Location','southwest', 'Box','off', 'FontSize',10)
    end
end

%% ============================================================
% Figure 2: Salinity profiles
%% ============================================================
fig2 = figure('Color','w','Units','pixels','Position',[80 70 1500 1050]);
tl2 = tiledlayout(fig2,3,3,'TileSpacing','compact','Padding','compact');

for ic = 1:nCast
    ax = nexttile(tl2,ic); hold(ax,'on'); box(ax,'on')

    plot(ax, obs_S_all{ic}, -obs_depth_all{ic}, 'k-', 'LineWidth', 2.3)

    for kk = 1:nSim
        plot(ax, mod_S_all{ic,kk}, -mod_depth_all{ic,kk}, '-', ...
            'Color', model_colors(kk,:), 'LineWidth', 1.5)
    end

    xlim(ax, salt_xlim)
    ylim(ax, [-depth_max_plot 0])
    set(ax,'FontSize',11,'LineWidth',1,'TickDir','out')

    title(ax, sprintf('Cast %d  %s', ic, datestr(cast_time(ic),'dd-mmm HH:MM')), ...
        'FontSize',12,'FontWeight','bold')

    if mod(ic-1,3) == 0
        ylabel(ax,'Depth (m)','FontSize',12)
    else
        set(ax,'YTickLabel',[])
    end

    if ic > 6
        xlabel(ax,'Salinity','FontSize',12)
    else
        set(ax,'XTickLabel',[])
    end

    if ic == 1
        legend(ax, [{'Obs'}, model_labels], 'Location','southwest', 'Box','off', 'FontSize',10)
    end
end

%% ============================================================
% Figure 3: Density profiles
%% ============================================================
fig3 = figure('Color','w','Units','pixels','Position',[110 80 1500 1050]);
tl3 = tiledlayout(fig3,3,3,'TileSpacing','compact','Padding','compact');

for ic = 1:nCast
    ax = nexttile(tl3,ic); hold(ax,'on'); box(ax,'on')

    plot(ax, obs_PD_all{ic}, -obs_depth_all{ic}, 'k-', 'LineWidth', 2.3)

    for kk = 1:nSim
        plot(ax, mod_PD_all{ic,kk}, -mod_depth_all{ic,kk}, '-', ...
            'Color', model_colors(kk,:), 'LineWidth', 1.5)
    end

    xlim(ax, pd_xlim)
    ylim(ax, [-depth_max_plot 0])
    set(ax,'FontSize',11,'LineWidth',1,'TickDir','out')

    title(ax, sprintf('Cast %d  %s', ic, datestr(cast_time(ic),'dd-mmm HH:MM')), ...
        'FontSize',12,'FontWeight','bold')

    if mod(ic-1,3) == 0
        ylabel(ax,'Depth (m)','FontSize',12)
    else
        set(ax,'YTickLabel',[])
    end

    if ic > 6
        xlabel(ax,'Density (kg m^{-3})','FontSize',12)
    else
        set(ax,'XTickLabel',[])
    end

    if ic == 1
        legend(ax, [{'Obs'}, model_labels], 'Location','southwest', 'Box','off', 'FontSize',10)
    end
end

%% ============================================================
% Per-profile RMSE tables
%% ============================================================
% profile_info = table( ...
%     cast_id, ...
%     datetime(cast_time,'ConvertFrom','datenum')', ...
%     cast_lat, cast_lon, ...
%     'VariableNames', {'Cast','Time','Lat','Lon'});

profile_info = table( ...
    cast_id(:), ...
    datetime(cast_time(:), 'ConvertFrom', 'datenum'), ...
    cast_lat(:), ...
    cast_lon(:), ...
    'VariableNames', {'Cast','Time','Lat','Lon'});

T_table = [profile_info, array2table(RMSE_T, 'VariableNames', ...
    {'GEN','SWOT_DA','Insitu_DA','Joint_DA'})];

S_table = [profile_info, array2table(RMSE_S, 'VariableNames', ...
    {'GEN','SWOT_DA','Insitu_DA','Joint_DA'})];

PD_table = [profile_info, array2table(RMSE_PD, 'VariableNames', ...
    {'GEN','SWOT_DA','Insitu_DA','Joint_DA'})];

disp(' ')
disp('Temperature RMSE by CTD cast (^oC):')
disp(T_table)

disp(' ')
disp('Salinity RMSE by CTD cast:')
disp(S_table)

disp(' ')
disp('Density RMSE by CTD cast (kg m^-3):')
disp(PD_table)

%% ============================================================
% Summary table: mean and std across 9 casts
%% ============================================================
summary_table = table( ...
    model_labels(:), ...
    mean(RMSE_T,  1, 'omitnan')', std(RMSE_T,  0, 1, 'omitnan')', ...
    mean(RMSE_S,  1, 'omitnan')', std(RMSE_S,  0, 1, 'omitnan')', ...
    mean(RMSE_PD, 1, 'omitnan')', std(RMSE_PD, 0, 1, 'omitnan')', ...
    'VariableNames', { ...
    'Simulation', ...
    'Temp_Mean', 'Temp_Std', ...
    'Salt_Mean', 'Salt_Std', ...
    'Density_Mean', 'Density_Std'});

disp(' ')
disp('Summary RMSE across 9 Leg 2 CTD casts:')
disp(summary_table)

%% ============================================================
% Optional save
%% ============================================================
% print(fig1, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/SI_CTD_Leg2_Temperature', '-dpng', '-r300');
% print(fig2, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/SI_CTD_Leg2_Salinity', '-dpng', '-r300');
% print(fig3, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/SI_CTD_Leg2_Density', '-dpng', '-r300');

%% ============================================================
% Local functions
%% ============================================================
function [lat_dec, lon_dec, t_num, t_str] = parse_ctd_header(cnv)

    lat_str = strtrim(cnv.instrumentheaders.NMEALatitude);
    lon_str = strtrim(cnv.instrumentheaders.NMEALongitude);
    t_str   = strtrim(cnv.instrumentheaders.NMEAUTCTime);

    tok = regexp(lat_str, '(\d+)\s+(\d+\.?\d*)\s*([NS])', 'tokens', 'once');
    lat_dec = str2double(tok{1}) + str2double(tok{2})/60;
    if strcmpi(tok{3}, 'S')
        lat_dec = -lat_dec;
    end

    tok = regexp(lon_str, '(\d+)\s+(\d+\.?\d*)\s*([EW])', 'tokens', 'once');
    lon_dec = str2double(tok{1}) + str2double(tok{2})/60;
    if strcmpi(tok{3}, 'W')
        lon_dec = -lon_dec;
    end

    t_utc = datetime(t_str, 'InputFormat', 'MMM dd yyyy HH:mm:ss', 'TimeZone', 'UTC');
    t_num = datenum(t_utc);
end

function depth_obs = get_ctd_depth(cnv)
    if isfield(cnv, 'depSM1')
        depth_obs = double(cnv.depSM1(:));
    elseif isfield(cnv, 'depSM')
        depth_obs = double(cnv.depSM(:));
    else
        error('No CTD depth field found.')
    end
end

function [best_name, best_it, best_time] = find_nearest_model_snapshot(model_path, dir_model, target_time)

    best_dt   = inf;
    best_name = '';
    best_it   = NaN;
    best_time = NaN;

    for ii = 1:length(dir_model)
        fname = fullfile(model_path, dir_model(ii).name);

        try
            ocean_time = ncread(fname, 'ocean_time');
        catch
            continue
        end

        tt = datenum(seconds(ocean_time) + datenum(1968,5,23));
        [dt_min, it] = min(abs(tt - target_time));

        if dt_min < best_dt
            best_dt   = dt_min;
            best_name = dir_model(ii).name;
            best_it   = it;
            best_time = tt(it);
        end
    end
end

function [depth_mod, T_prof, S_prof, PD_prof] = extract_model_profile(model, it_model, lon_cast, lat_cast)

    depth_mod = double(model.depth(:));
    nZ = length(depth_mod);

    T_prof  = nan(nZ,1);
    S_prof  = nan(nZ,1);
    PD_prof = nan(nZ,1);

    for iz = 1:nZ
        T_prof(iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.temp(:,:,iz,it_model))', ...
            lon_cast, lat_cast, 'linear');

        S_prof(iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.salt(:,:,iz,it_model))', ...
            lon_cast, lat_cast, 'linear');

        PD_prof(iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.rho(:,:,iz,it_model))', ...
            lon_cast, lat_cast, 'linear');
    end
end
