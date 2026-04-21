% ============================================================
% Figure 3 release script
% Panels: Figure 3a-o
% Purpose: Across-eddy SSH, temperature, and salinity sections along the glider transect.
%
% Notes for release:
% - Run this script directly.
% - Key user choices are section_id, file_id, and time_id near the top.
% - This script generates the glider/model section panels and RMSE values used in Figure 3.
% - Keep figure settings unchanged unless you are intentionally updating the published figure.
% ============================================================


clear all; close all; clc

%% ============================================================
% Toolboxes and colormaps
% ============================================================
addpath TEOS_10
addpath TEOS_10/library/
addpath Sea-Bird-Toolbox-master/CTD_CNV/
addpath cmocean

RT = load("rt_colormaps.mat");

%% ============================================================
% Paths
% ============================================================
glider_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/';
SWOT_path   = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';

WMOP_paths = { ...
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v6/', ... % GEN
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7/', ... % SWOT DA
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8/', ... % In-situ DA
    '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9/'  ... % Joint DA
    };

version_labels = {'GEN','SWOT DA','In-situ DA','Joint DA'};

colors =[ 0.0660    0.4430    0.7450
          0.2310    0.6660    0.1960
          0.5210    0.0860    0.8190
          0.9608    0.4667    0.1608 ];

set(groot,'defaultAxesFontName','Helvetica')
set(groot,'defaultTextFontName','Helvetica')

%% ============================================================
% Load glider data
% ============================================================
[~, glider2] = read_nc_file_struct([glider_path 'dep0007_sdeep09_scb-sldeep009_L2_2023-04-25_data_dt.nc']);
glider2.Time = datenum(seconds(glider2.time) + datenum(1970,01,01));

glider2.CT = gsw_CT_from_t(glider2.salinity, glider2.temperature, 0);
glider2.SA = gsw_SA_from_SP(glider2.salinity, 0, glider2.longitude, glider2.latitude);
glider2.PD = gsw_rho(glider2.SA, glider2.CT, 0) - 1000;

%% ============================================================
% Define glider sections
% ============================================================
g2{1} = find(glider2.Time >= datenum(2023,4,28,9,48,00) ...
          & glider2.Time <= datenum(2023,5,1,16,50,00));

g2{2} = find(glider2.Time >= datenum(2023,5,1,16,50,00) ...
          & glider2.Time <= datenum(2023,5,3,15,45,00));

% g2{3} = find(glider2.Time >= datenum(2023,5,3,15,45,00) ...
%           & glider2.Time <= datenum(2023,5,7,07,04,00));

g2{3} = find(glider2.Time >= datenum(2023,5,3,15,45,00) ...
          & glider2.Time <= datenum(2023,5,7,00,00,00));

glider_hr = [2, 6, 6];

%% ============================================================
% User choices
% ============================================================
section_id = 3;      % use g2{3}
file_id    = 16;     % model file index
time_id    = glider_hr(section_id);
ind_sec    = g2{section_id};

lat_tran   = glider2.latitude(ind_sec);
lon_tran   = glider2.longitude(ind_sec);

lat_xlim   = [39.54 40.15];
depth_ylim = [-600 0];

temp_lim   = [13.1 13.8];
salt_lim   = [38.2 38.6];
ssh_lim    = [-2 2];

rho_cont_all   = 28.8:0.05:29.1;
rho_cont_label = [28.9 28.95];

panel_letters = arrayfun(@(k) ['(' char('a'+k-1) ')'], 1:15, 'UniformOutput', false);

%% ============================================================
% Load SWOT data and interpolate to glider track
% ============================================================
dir_SWOT = dir(fullfile(SWOT_path, '*.nc'));
[~, data_swot] = read_nc_file_struct(fullfile(SWOT_path, dir_SWOT(34).name));
data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));

valid_swot = isfinite(data_swot.longitude) & isfinite(data_swot.latitude) & isfinite(data_swot.ssha_filtered);

Fswot = scatteredInterpolant( ...
    data_swot.longitude(valid_swot), ...
    data_swot.latitude(valid_swot), ...
    data_swot.ssha_filtered(valid_swot), ...
    'natural', 'none');

swot_ssh_line = Fswot(lon_tran, lat_tran);

%% ============================================================
% Containers
% ============================================================
SSH_model = cell(1,4);
T_model   = cell(1,4);
S_model   = cell(1,4);
PD_model  = cell(1,4);
RMSE_ssh  = nan(1,4);

%% ============================================================
% Read model data and interpolate onto glider transect
% ============================================================
for kk = 1:4

    dir_WMOP = dir(fullfile(WMOP_paths{kk}, '*.nc'));
    [~, model] = read_nc_file_struct(fullfile(WMOP_paths{kk}, dir_WMOP(file_id).name));

    model.time = datenum(seconds(model.ocean_time) + datenum(1968,5,23));
    model.rho  = gsw_rho(model.salt, model.temp, 40) - 1000;

    if kk == 1
        model_depth = double(model.depth(:));
        model_time  = model.time(time_id);
    end

    nZ = length(model_depth);
    nP = length(ind_sec);

    T_sec  = nan(nP, nZ);
    S_sec  = nan(nP, nZ);
    PD_sec = nan(nP, nZ);

    for iz = 1:nZ
        PD_sec(:,iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.rho(:,:,iz,time_id))', lon_tran, lat_tran, 'linear');

        T_sec(:,iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.temp(:,:,iz,time_id))', lon_tran, lat_tran, 'linear');

        S_sec(:,iz) = interp2(double(model.lon), double(model.lat), ...
            squeeze(model.salt(:,:,iz,time_id))', lon_tran, lat_tran, 'linear');
    end

    PD_model{kk} = PD_sec;
    T_model{kk}  = T_sec;
    S_model{kk}  = S_sec;

    zeta = squeeze(model.zeta(:,:,time_id))';

    axis_def = [0.8 3 39 41];
    lon_min = axis_def(1);
    lon_max = axis_def(2);
    lat_min = axis_def(3);
    lat_max = axis_def(4);

    [Lon, Lat] = meshgrid(model.lon, model.lat);
    mask = Lon >= lon_min & Lon <= lon_max & Lat >= lat_min & Lat <= lat_max;

    zeta_masked = zeta;
    zeta_masked(~mask) = NaN;
    zeta_anom = zeta - nanmean(zeta_masked(:));

    ssh_line = interp2(double(model.lon), double(model.lat), zeta_anom, lon_tran, lat_tran, 'linear');
    SSH_model{kk} = ssh_line;

    A = (ssh_line      - nanmean(ssh_line))*100;
    B = (swot_ssh_line - nanmean(swot_ssh_line))*100;
    RMSE_ssh(kk) = sqrt(nanmean((A - B).^2));
end

%% ============================================================
% Glider data on model depth grid
% only used for SSH and any comparison if needed
%% ============================================================
nP = length(ind_sec);
nZ = length(model_depth);

PD_glider_on_model = nan(nP, nZ);
T_glider_on_model  = nan(nP, nZ);
S_glider_on_model  = nan(nP, nZ);

for ip = 1:nP
    PD_glider_on_model(ip,:) = interp1(glider2.depth, glider2.PD(:,ind_sec(ip)), ...
        model_depth, 'linear', NaN);

    T_glider_on_model(ip,:)  = interp1(glider2.depth, glider2.temperature(:,ind_sec(ip)), ...
        model_depth, 'linear', NaN);

    S_glider_on_model(ip,:)  = interp1(glider2.depth, glider2.salinity(:,ind_sec(ip)), ...
        model_depth, 'linear', NaN);
end

%% ============================================================
% RMSE of model T/S against glider on common grid
% place this right after the block where
% T_glider_on_model and S_glider_on_model are computed
%% ============================================================

RMSE_T = nan(1,4);
RMSE_S = nan(1,4);

RMSE_T = nan(1,4);
RMSE_S = nan(1,4);

iz_rmse = find(model_depth >= 0 & model_depth <= 600);

for kk = 1:4
    A_T = T_model{kk}(:, iz_rmse);
    B_T = T_glider_on_model(:, iz_rmse);

    A_S = S_model{kk}(:, iz_rmse);
    B_S = S_glider_on_model(:, iz_rmse);

    RMSE_T(kk) = sqrt(nanmean((A_T(:) - B_T(:)).^2));
    RMSE_S(kk) = sqrt(nanmean((A_S(:) - B_S(:)).^2));
end


%% ============================================================
% Real glider data on native glider depth grid
% these will be used in panels (j) and (o)
%% ============================================================
[lat_unique, idx_unique] = unique(lat_tran, 'stable');

lat_glider_real = lat_unique;
ind_unique      = ind_sec(idx_unique);

T_glider_real  = glider2.temperature(:, ind_unique);
S_glider_real  = glider2.salinity(:, ind_unique);
PD_glider_real = glider2.PD(:, ind_unique);


swot_ssh_u = swot_ssh_line(idx_unique);

%% ============================================================
% Dynamic height SSH from glider
% ============================================================
p = gsw_p_from_z(-glider2.depth, nanmean(lat_tran));
p_ref = 700;

dyn_geo = gsw_geo_strf_dyn_height(glider2.SA(:,ind_sec), glider2.CT(:,ind_sec), p, p_ref);
dyn_h   = dyn_geo ./ 9.81;

ssh_glider = dyn_h(16,:);
ssh_glider = fillmissing(ssh_glider, 'linear');
ssh_glider = smoothdata(ssh_glider, 'gaussian', 4);
ssh_glider_u = ssh_glider(idx_unique);

%% ============================================================
% Model time marker on glider transect
%% ============================================================
[~, idx_tstar] = min(abs(glider2.Time(ind_sec) - model_time));
lat_star = lat_tran(idx_tstar);

%% ============================================================
% Colormap
% ============================================================
% if isfield(RT,'rt_colormaps') && isfield(RT.rt_colormaps,'section')
%     cmap_section = RT.rt_colormaps.section;
% % else
% %     cmap_section = parula(256);
% end

%% ============================================================
% Manual layout for unequal row heights
% Top row is shorter
%% ============================================================
fig = figure('Color','w','Units','pixels','Position',[40 40 2050 1150]);

margin_left   = 0.055;
margin_right  = 0.09;
margin_top    = 0.100;
margin_bottom = 0.10;

gap_x = 0.018;
gap_y = 0.026;

row_h_top = 0.14;
row_h_mid = 0.33* 0.85;
row_h_bot = 0.33* 0.85;

col_w = (1 - margin_left - margin_right - 4*gap_x) / 5;

% y_bot = margin_bottom;
% y_mid = y_bot + row_h_bot + gap_y;
% y_top = y_mid + row_h_mid + gap_y;

% redistribute vertically so the whole figure stays balanced
total_h_used = row_h_top + row_h_mid + row_h_bot + 2*gap_y;
free_h = 1 - margin_top - margin_bottom - total_h_used;

% put half of the freed space between top and middle,
% and half between middle and bottom
extra_gap = free_h / 2;

y_bot = margin_bottom;
y_mid = y_bot + row_h_bot + gap_y + extra_gap;
y_top = y_mid + row_h_mid + gap_y + extra_gap;

% manual vertical tweaks
dy_top = -0.015;   % move top row downward
dy_bot = +0.015;   % move bottom row upward

y_top = y_top + dy_top;
y_bot = y_bot + dy_bot;

xpos = zeros(1,5);
for cc = 1:5
    xpos(cc) = margin_left + (cc-1)*(col_w + gap_x);
end

% sgtitle(sprintf('Across-eddy section, %s', datestr(model_time,'dd-mmm-yyyy HH:MM')), ...
%     'FontSize',16,'FontWeight','bold')
annotation('textbox', [0.32 0.915 0.36 0.035], ...
    'String', sprintf('Across-eddy section, %s', datestr(model_time,'dd-mmm-yyyy HH:MM')), ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

% cmap_temp = cmocean('thermal',256);
cmap_temp = RT.rt_colormaps.greenbrown; 
cmap_salt = RT.rt_colormaps.sst2;
% cmap_salt = RT.rt_colormaps.greenbrown;
%% ============================================================
% TOP ROW: SSH (a-e)
%% ============================================================
for kk = 1:4
    ax = axes('Position',[xpos(kk) y_top col_w row_h_top]); hold(ax,'on'); box(ax,'on')

    plot(ax, lat_tran, (SSH_model{kk} - nanmean(SSH_model{kk}))*100, ...
        'Color', colors(kk,:), 'LineWidth', 2.8)

    plot(ax, lat_tran, (swot_ssh_line - nanmean(swot_ssh_line))*100, ...
        'k', 'LineWidth', 2.4)

    xlim(ax, lat_xlim)
    ylim(ax, ssh_lim)
    grid(ax,'on')
    set(ax,'FontSize',12,'LineWidth',1,'TickDir','out','XTickLabel',[])

    title(ax, version_labels{kk}, 'FontSize',14, 'FontWeight','bold')

    if kk == 1
        ylabel(ax,'SSH (cm)','FontSize',14)
        lgd = legend(ax, {'Model','SWOT'}, 'Location','northeast');
        lgd.Box = 'off';
    else
        set(ax,'YTickLabel',[])
    end

    text(ax, 0.02, 0.94, panel_letters{kk}, 'Units','normalized', ...
        'FontSize',14, 'FontWeight','bold')

    text(ax, 0.30, 0.08, sprintf('RMSE = %.2f cm', RMSE_ssh(kk)), ...
        'Units','normalized', 'FontSize',14, 'FontWeight','bold', ...
        'Color', colors(kk,:))
end

% Panel (e): Glider SSH vs SWOT SSH
ax_e = axes('Position',[xpos(5) y_top col_w row_h_top]); hold(ax_e,'on'); box(ax_e,'on')

plot(ax_e, lat_glider_real, (ssh_glider_u - nanmean(ssh_glider_u))*100, ...
    'b', 'LineWidth', 2.8)

plot(ax_e, lat_glider_real, (swot_ssh_u - nanmean(swot_ssh_u))*100, ...
    'k', 'LineWidth', 2.4)

% plot(ax_e, lat_star, ssh_lim(2)-0.12, 'k*', 'MarkerSize', 9, 'LineWidth', 1.2)

xlim(ax_e, lat_xlim)
ylim(ax_e, ssh_lim)
grid(ax_e,'on')
set(ax_e,'FontSize',12,'LineWidth',1,'TickDir','out','XTickLabel',[],'YTickLabel',[])

glider_time_str = [datestr(glider2.Time(ind_sec(1)), 'dd-mmm HH:MM'), ...
                   ' to ', ...
                   datestr(glider2.Time(ind_sec(end)), 'dd-mmm HH:MM')];


% title(ax_e, ['Glider ' glider_time_str], 'FontSize',14, 'FontWeight','bold')

title(ax_e, ['Glider 3 May 17:09 -- 6 May 23:38'], 'FontSize',14, 'FontWeight','bold')


% lgd = legend(ax_e, {'Glider','SWOT'}, 'Location','northeast');
% lgd.Box = 'off';

text(ax_e, 0.02, 0.94, panel_letters{5}, 'Units','normalized', ...
    'FontSize',14, 'FontWeight','bold')

%% ============================================================
% MIDDLE ROW: Temperature (f-j)
%% ============================================================
for kk = 1:4
    ax = axes('Position',[xpos(kk) y_mid col_w row_h_mid]); hold(ax,'on'); box(ax,'on')

    pcolor(ax, lat_tran, -model_depth, T_model{kk}'); shading(ax,'interp')
    

    colormap(ax, cmap_temp)
    caxis(ax, temp_lim)

    contour(ax, lat_tran, -model_depth, PD_model{kk}', rho_cont_all, ...
        'k', 'LineWidth', 0.45)

    [Ctmp, htmp] = contour(ax, lat_tran, -model_depth, PD_model{kk}', rho_cont_label, ...
        'k', 'LineWidth', 0.9);
    clabel(Ctmp, htmp, 'FontSize', 12, 'Color', 'k', 'LabelSpacing', 350)

    xlim(ax, lat_xlim)
    ylim(ax, depth_ylim)
    set(ax,'FontSize',12,'LineWidth',1,'TickDir','out','XTickLabel',[])

    if kk == 1
        ylabel(ax,'Depth (m)','FontSize',14)
    else
        set(ax,'YTickLabel',[])
    end

    text(ax, 0.02, 0.96, panel_letters{5+kk}, 'Units','normalized', ...
        'FontSize',14, 'FontWeight','bold', 'Color','w')

    text(ax, 0.3, 0.06, sprintf('RMSE = %.2f ^oC', RMSE_T(kk)), ...
    'Units','normalized', 'FontSize',14, 'FontWeight','bold', ...
    'Color','k', 'Margin',2)

end

% Panel (j): real glider temperature on native depth grid
ax_j = axes('Position',[xpos(5) y_mid col_w row_h_mid]); hold(ax_j,'on'); box(ax_j,'on')

pcolor(ax_j, lat_glider_real, -glider2.depth, T_glider_real); shading(ax_j,'flat')
colormap(ax_j, cmap_temp)
caxis(ax_j, temp_lim)

contour(ax_j, lat_glider_real, -glider2.depth, PD_glider_real, rho_cont_all, ...
    'k', 'LineWidth', 0.45)

[Cj, hj] = contour(ax_j, lat_glider_real, -glider2.depth, PD_glider_real, rho_cont_label, ...
    'k', 'LineWidth', 0.9);
clabel(Cj, hj, 'FontSize', 12, 'Color', 'k', 'LabelSpacing', 350)

% plot(ax_j, lat_star, -8, 'k*', 'MarkerSize', 9, 'LineWidth', 1.2)

xlim(ax_j, lat_xlim)
ylim(ax_j, depth_ylim)
set(ax_j,'FontSize',12,'LineWidth',1,'TickDir','out','XTickLabel',[],'YTickLabel',[])

text(ax_j, 0.02, 0.96, panel_letters{10}, 'Units','normalized', ...
    'FontSize',14, 'FontWeight','bold', 'Color','w')

cb1 = colorbar(ax_j, 'eastoutside');
cb1.Position = [xpos(5)+col_w+0.005, y_mid, 0.010, row_h_mid];
cb1.Label.String = 'Temperature (°C)';
cb1.FontSize = 11;
cb1.TickDirection = 'out';
cb1.YAxisLocation = 'Left';

%% ============================================================
% BOTTOM ROW: Salinity (k-o)
%% ============================================================
for kk = 1:4
    ax = axes('Position',[xpos(kk) y_bot col_w row_h_bot]); hold(ax,'on'); box(ax,'on')

    pcolor(ax, lat_tran, -model_depth, S_model{kk}'); shading(ax,'interp')
    % colormap(ax, cmap_section)
    % colormap(ax,RT.rt_colormaps.sst)
    colormap(ax, cmap_salt)
    caxis(ax, salt_lim)

    contour(ax, lat_tran, -model_depth, PD_model{kk}', rho_cont_all, ...
        'k', 'LineWidth', 0.45)

    [Csal, hsal] = contour(ax, lat_tran, -model_depth, PD_model{kk}', rho_cont_label, ...
        'k', 'LineWidth', 0.9);
    clabel(Csal, hsal, 'FontSize', 12, 'Color', 'k', 'LabelSpacing', 350)

    xlim(ax, lat_xlim)
    ylim(ax, depth_ylim)
    set(ax,'FontSize',12,'LineWidth',1,'TickDir','out')

    if kk == 1
        ylabel(ax,'Depth (m)','FontSize',14)
    else
        set(ax,'YTickLabel',[])
    end

    xlabel(ax,'Lat (\circN)','FontSize',14)

    text(ax, 0.02, 0.96, panel_letters{10+kk}, 'Units','normalized', ...
        'FontSize',14, 'FontWeight','bold', 'Color','w')

    text(ax, 0.3, 0.06, sprintf('RMSE = %.3f g kg^{-1}', RMSE_S(kk)), ...
        'Units','normalized', 'FontSize',14, 'FontWeight','bold', ...
        'Color','w','Margin',2)
end

% Panel (o): real glider salinity on native depth grid
ax_o = axes('Position',[xpos(5) y_bot col_w row_h_bot]); hold(ax_o,'on'); box(ax_o,'on')

pcolor(ax_o, lat_glider_real, -glider2.depth, S_glider_real); shading(ax_o,'flat')
colormap(ax_o, cmap_salt)
% colormap(ax_o,RT.rt_colormaps.sst)

% colormap(ax_o,RT.rt_colormaps.greenbrown); caxis(ax_o, [38.2 38.6])
% cmocean('curl'); caxis(ax_o, [38.2 38.6])


caxis(ax_o, salt_lim)

contour(ax_o, lat_glider_real, -glider2.depth, PD_glider_real, rho_cont_all, ...
    'k', 'LineWidth', 0.45)

[Co, ho] = contour(ax_o, lat_glider_real, -glider2.depth, PD_glider_real, rho_cont_label, ...
    'k', 'LineWidth', 0.9);
clabel(Co, ho, 'FontSize', 12, 'Color', 'k', 'LabelSpacing', 350)

plot(ax_o, lat_star, -8, 'k*', 'MarkerSize', 9, 'LineWidth', 1.2)

xlim(ax_o, lat_xlim)
ylim(ax_o, depth_ylim)
set(ax_o,'FontSize',12,'LineWidth',1,'TickDir','out','YTickLabel',[])

xlabel(ax_o,'Lat (°N)','FontSize',14)

text(ax_o, 0.02, 0.96, panel_letters{15}, 'Units','normalized', ...
    'FontSize',14, 'FontWeight','bold', 'Color','w')

cb2 = colorbar(ax_o, 'eastoutside');
cb2.Position = [xpos(5)+col_w+0.005, y_bot, 0.010, row_h_bot];
cb2.YAxisLocation = 'Left';
cb2.Label.String = 'Salinity (g kg^{-1})';
cb2.FontSize = 11;
cb2.TickDirection = 'out';

%% ============================================================
print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/GRL_Fig3'], '-dpng', '-r300');      
