% ============================================================
% Figure 4 release script
% Panels: Figure 4a-c
% Purpose: 2D SSH maps on 9 May for SWOT DA, Joint DA, and SWOT, with ADCP tracks.
%
% Notes for release:
% - Run this script first for Figure 4 map panels.
% - This script corresponds only to panels (a-c).
% - Figure 4 was assembled from multiple scripts, listed in the workflow note.
% ============================================================


clear; close all; clc

%% ============================================================
% 1 x 3 figure on 09 May 2023
% (a) SWOT DA
% (b) Joint DA
% (c) SWOT
% Keep ADCP tracks on all panels
%% ============================================================

addpath cmocean
addpath TEOS_10
addpath TEOS_10/library/

RT = load("rt_colormaps.mat");
C  = load('coastline_full_westmed_nolakes.mat');

if isfield(C,'ncst')
    ncst = C.ncst;
else
    fn = fieldnames(C);
    ncst = C.(fn{1});
end

colors =[ 0.0660    0.4430    0.7450
          0.2310    0.6660    0.1960
          0.5210    0.0860    0.8190
          0.9608    0.4667    0.1608          
          ];


%% ============================================================
% Paths
%% ============================================================
WMOP_path_swotda = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7_fr_07_05/'; % SWOT DA
WMOP_path_joint  = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9_fr_07_05/'; % Joint DA
SWOT_path        = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';
adcp_path        = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/ADCP/';   % revise if needed

dir_swotda = dir(fullfile(WMOP_path_swotda,'*.nc'));
dir_joint  = dir(fullfile(WMOP_path_joint ,'*.nc'));
dir_swot   = dir(fullfile(SWOT_path,'*.nc'));

%% ============================================================
% Load ADCP data and define target transects
%% ============================================================
[~, data_Leg] = read_nc_file_struct([adcp_path 'ADCP_OS150_FaSt-SWOT_LEG2.nc']);
data_Leg.Time = data_Leg.time + datenum(2023,1,1);

ind{1} = find(data_Leg.Time > datenum(2023,05,08,20,30,00) & ...
              data_Leg.Time < datenum(2023,05,09,03,35,00));

ind{2} = find(data_Leg.Time > datenum(2023,05,09,04,40,00) & ...
              data_Leg.Time < datenum(2023,05,09,14,30,00));

adcp1_lon = data_Leg.lon(ind{1});
adcp1_lat = data_Leg.lat(ind{1});

adcp2_lon = data_Leg.lon(ind{2});
adcp2_lat = data_Leg.lat(ind{2});




CTD_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/CTDs/Leg2/';

CTD_dir = dir(CTD_path);
% ustation: upcast, dstation: downcast

k = 1;
Z(1:3150,1:(12 -3)) = NaN;
T(1:3150,1:(12 -3)) = NaN;
T_CTD(1:3150,1:(12 -3)) = NaN;


for i = 4 : 12 %chosen downcast files only, length(CTD_dir)
 
        cnv = readSBScnv([CTD_path CTD_dir(i).name]);
        % --- grab raw strings from the struct
        lat_str = cnv.instrumentheaders.NMEALatitude;   % e.g. '39 41.21 N'
        lon_str = cnv.instrumentheaders.NMEALongitude;  % e.g. '001 26.46 E'
        t_str   = cnv.instrumentheaders.NMEAUTCTime;    % e.g. 'Apr 26 2023 10:58:09'
        
%% check that SA make comparison worse
        Z(1:length(cnv.depSM1),k) = cnv.depSM1;                         % depth vector
        P = gsw_p_from_z(-cnv.depSM1,40);                           
        % S(1:length(cnv.t090C),k)  = gsw_SA_from_SP(cnv.psal0,P,1.55,39.7);
        S  = cnv.psal0; %gsw_SA_from_SP(cnv.psal0,P,1.55,39.7);
        % T(1:length(cnv.t090C),k)  = gsw_CT_from_t(S,cnv.t090C,P);  
        T = cnv.t090C; 
        % rho(1:length(cnv.t090C),k)  = gsw_rho(S,T,0);        
        % T = rho;        
        T_CTD(1:length(T),i-3) = T;



        % --- grab raw strings from the struct     
        t_str   = cnv.instrumentheaders.NMEAUTCTime;    % e.g. 'Apr 26 2023 10:58:09'
        CTD_time{i-3} = t_str;

        % --- convert 'dd mm.mmm H' / 'ddd mm.mmm H' -> decimal degrees
        tok = regexp(lat_str,'(\d+)\s+(\d+\.?\d*)\s*([NS])','tokens','once');
        CTD_lat(k) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='S', lat = -lat; end
        
        tok = regexp(lon_str,'(\d+)\s+(\d+\.?\d*)\s*([EW])','tokens','once');
        CTD_lon(k) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='W', lon = -lon; end
       
        % --- UTC time as datetime (and datenum if you need it)
        t_utc = datetime(t_str,'InputFormat','MMM dd yyyy HH:mm:ss','TimeZone','UTC');
        time_all(k) = datenum(t_utc);   % optional, legacy

        k = k + 1 ;
    end

glider_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/';

% Load glider data
%% ============================================================
[~, glider1] = read_nc_file_struct([glider_path 'dep0040_sdeep01_scb-sldeep001_L2_2023-04-25_data_dt.nc']);
glider1.Time = datenum(seconds(glider1.time) + datenum(1970,01,01));
% glider1.CT = gsw_CT_from_t(glider1.salinity, glider1.temperature, 0);
% glider1.SA = gsw_SA_from_SP(glider1.salinity, 0, glider1.longitude, glider1.latitude);
% glider1.PD = gsw_rho(glider1.SA, glider1.CT, 0) - 1000;

[~, glider2] = read_nc_file_struct([glider_path 'dep0007_sdeep09_scb-sldeep009_L2_2023-04-25_data_dt.nc']);
glider2.Time = datenum(seconds(glider2.time) + datenum(1970,01,01));
% glider2.CT = gsw_CT_from_t(glider2.salinity, glider2.temperature, 0);
% glider2.SA = gsw_SA_from_SP(glider2.salinity, 0, glider2.longitude, glider2.latitude);
% glider2.PD = gsw_rho(glider2.SA, glider2.CT, 0) - 1000;

%% ============================================================
% Define glider sections
%% ============================================================
g1_ind = find(glider1.Time >= datenum(2023,5,7,23,40,00) & ...
              glider1.Time <= datenum(2023,5,10,07,03,00));

g2_ind = find(glider2.Time >= datenum(2023,5,7,07,38,00) & ...
              glider2.Time <= datenum(2023,5,10,07,12,00));

%% ============================================================
% Settings
%% ============================================================
target_date    = datetime(2023,5,9);

axis_def       = [0.8 2.2 39.2 40.7];
axis_def_mask  = [0.8 3.0 39.0 41.0];

lon_target     = 2.0;
lat_target     = 40.0;

cax_lim        = [-7 7];
cont_lev       = -10:1:10;

lon_min = axis_def_mask(1);
lon_max = axis_def_mask(2);
lat_min = axis_def_mask(3);
lat_max = axis_def_mask(4);

% panel_letters = {'(a)','(b)','(c)'};

%% ============================================================
% Precompute SWOT representative time and day
% Representative time = closest point in swath to (2E, 40N)
%% ============================================================
swot_closest_time = nan(length(dir_swot),1);
swot_day          = nan(length(dir_swot),1);

for k = 1:length(dir_swot)
    fname = fullfile(SWOT_path, dir_swot(k).name);

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
% Pick SWOT file for 09 May
%% ============================================================
target_dnum = datenum(target_date);
idx_swot_day = find(swot_day == floor(target_dnum));

if isempty(idx_swot_day)
    valid_swot = find(isfinite(swot_closest_time));
    [~, jj] = min(abs(swot_closest_time(valid_swot) - target_dnum));
    idx_swot = valid_swot(jj);
else
    [~, jj] = min(abs(swot_closest_time(idx_swot_day) - target_dnum));
    idx_swot = idx_swot_day(jj);
end

swot_fname = fullfile(SWOT_path, dir_swot(idx_swot).name);
[~, data_swot] = read_nc_file_struct(swot_fname);
data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));

t_swot = swot_closest_time(idx_swot);

%% ============================================================
% SWOT masked SSH anomaly
%% ============================================================
mask_swot = data_swot.longitude >= lon_min & data_swot.longitude <= lon_max & ...
            data_swot.latitude  >= lat_min & data_swot.latitude  <= lat_max;

ssha_masked = data_swot.ssha_filtered;
ssha_masked(~mask_swot) = NaN;

swot_ssh = (ssha_masked - nanmean(ssha_masked(:))) * 100;

% same cropping as your old script
irow0 = 36;

%% ============================================================
% Precompute model file day metadata
%% ============================================================
swotda_file_day = nan(length(dir_swotda),1);
joint_file_day  = nan(length(dir_joint),1);

for k = 1:length(dir_swotda)
    fname = fullfile(WMOP_path_swotda, dir_swotda(k).name);
    try
        ocean_time = ncread(fname,'ocean_time');
        tt = datenum(seconds(ocean_time) + datenum(1968,5,23));
        swotda_file_day(k) = floor(tt(1));
    catch
        swotda_file_day(k) = NaN;
    end
end

for k = 1:length(dir_joint)
    fname = fullfile(WMOP_path_joint, dir_joint(k).name);
    try
        ocean_time = ncread(fname,'ocean_time');
        tt = datenum(seconds(ocean_time) + datenum(1968,5,23));
        joint_file_day(k) = floor(tt(1));
    catch
        joint_file_day(k) = NaN;
    end
end

%% ============================================================
% Read SWOT DA model at closest time to SWOT
%% ============================================================
idx_swotda_day = find(swotda_file_day == floor(target_dnum));

if isempty(idx_swotda_day)
    valid_swotda = find(isfinite(swotda_file_day));
    [~, jj] = min(abs(swotda_file_day(valid_swotda) - floor(target_dnum)));
    idx_swotda = valid_swotda(jj);
else
    idx_swotda = idx_swotda_day(1);
end

swotda_fname = fullfile(WMOP_path_swotda, dir_swotda(idx_swotda).name);
[~, data_swotda] = read_nc_file_struct(swotda_fname);
data_swotda.time = datenum(seconds(data_swotda.ocean_time) + datenum(1968,5,23));

[~, it_swotda] = min(abs(data_swotda.time - t_swot));
t_swotda = data_swotda.time(it_swotda);

zeta_swotda = data_swotda.zeta(:,:,it_swotda)';

[Lon_swotda, Lat_swotda] = meshgrid(data_swotda.lon, data_swotda.lat);
mask_swotda = Lon_swotda >= lon_min & Lon_swotda <= lon_max & ...
              Lat_swotda >= lat_min & Lat_swotda <= lat_max;

zeta_swotda_masked = zeta_swotda;
zeta_swotda_masked(~mask_swotda) = NaN;

zeta_swotda_anom = zeta_swotda - nanmean(zeta_swotda_masked(:));

%% ============================================================
% Read Joint DA model at closest time to SWOT
%% ============================================================
idx_joint_day = find(joint_file_day == floor(target_dnum));

if isempty(idx_joint_day)
    valid_joint = find(isfinite(joint_file_day));
    [~, jj] = min(abs(joint_file_day(valid_joint) - floor(target_dnum)));
    idx_joint = valid_joint(jj);
else
    idx_joint = idx_joint_day(1);
end

joint_fname = fullfile(WMOP_path_joint, dir_joint(idx_joint).name);
[~, data_joint] = read_nc_file_struct(joint_fname);
data_joint.time = datenum(seconds(data_joint.ocean_time) + datenum(1968,5,23));

[~, it_joint] = min(abs(data_joint.time - t_swot));
t_joint = data_joint.time(it_joint);

zeta_joint = data_joint.zeta(:,:,it_joint)';

[Lon_joint, Lat_joint] = meshgrid(data_joint.lon, data_joint.lat);
mask_joint = Lon_joint >= lon_min & Lon_joint <= lon_max & ...
             Lat_joint >= lat_min & Lat_joint <= lat_max;

zeta_joint_masked = zeta_joint;
zeta_joint_masked(~mask_joint) = NaN;

zeta_joint_anom = zeta_joint - nanmean(zeta_joint_masked(:));

%% ============================================================
% SWOT swath polygon for model panels
%% ============================================================
L1_lon = data_swot.longitude(40,:);
L1_lat = data_swot.latitude(40,:);

L2_lon = data_swot.longitude(end-4,:);
L2_lat = data_swot.latitude(end-4,:);

poly_lon = [L1_lon fliplr(L2_lon) L1_lon(1)];
poly_lat = [L1_lat fliplr(L2_lat) L1_lat(1)];

%% ============================================================
% Figure
%% ============================================================
fig = figure('Color','w','Units','pixels','Position',[90 120 1150 460]);
tl  = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    % fig = figure('Color','w','Units','pixels','Position',[100 120 1250 350]);
    % tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

colormap(fig, RT.rt_colormaps.section);

%% ------------------------------------------------------------
% (a) SWOT DA
%% ------------------------------------------------------------
ax1 = nexttile(tl,1); hold(ax1,'on')

hp = pcolor(ax1, data_swotda.lon, data_swotda.lat, zeta_swotda_anom*100);
shading(ax1,'flat')

[c1,b1] = contour(ax1, data_swotda.lon, data_swotda.lat, zeta_swotda_anom*100, ...
    cont_lev, 'k', 'LineWidth', 0.8);         
clabel(c1,b1,'fontsize',14)

plot(ax1, ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)
plot(ax1, poly_lon, poly_lat, '-', 'Color', [1 1 1], 'LineWidth', 2.0)

% ===== SWOT TRACKS USED FOR MASK =====
L1_lon = data_swot.longitude(39+1,:);
L1_lat = data_swot.latitude(39+1,:);
L2_lon = data_swot.longitude(end-4,:);
L2_lat = data_swot.latitude(end-4,:);

% ===== BUILD CLOSED POLYGON =====
poly_lon = [L1_lon fliplr(L2_lon) L1_lon(1)];
poly_lat = [L1_lat fliplr(L2_lat) L1_lat(1)];

% ===== PLOT CLOSED POLYGON =====
plot(ax1,poly_lon, poly_lat, '-', 'Color', colors(2,:), 'LineWidth', 2.2)

% ===== CREATE MASK =====
[Lon2D, Lat2D] = meshgrid(data_swotda.lon, data_swotda.lat);
in = inpolygon(Lon2D, Lat2D, poly_lon, poly_lat);

alpha_mask = 0.35 * ones(size(in));
alpha_mask(in) = 1;

% ===== APPLY TRANSPARENCY =====
set(hp,'AlphaData',alpha_mask);
set(hp,'FaceAlpha','flat');
set(hp,'AlphaDataMapping','none');


% ADCP transect 2: solid
% plot(ax1, adcp2_lon, adcp2_lat, 'w-', 'LineWidth', 3.0)
plot(ax1, adcp2_lon, adcp2_lat, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2)
% CTD stations
plot(ax1, CTD_lon, CTD_lat,'o','MarkerFaceColor',[0.5210    0.0860    0.8190],'MarkerSize',6)
% Glider track
plot(ax1, glider1.longitude(g1_ind), glider1.latitude(g1_ind),'b-', 'LineWidth', 2)
plot(ax1, glider2.longitude(g2_ind), glider2.latitude(g2_ind),'b-', 'LineWidth', 2)

axis(ax1,'equal')
axis(ax1, axis_def)
box(ax1,'on')
caxis(ax1, cax_lim)

set(ax1,'FontSize',12,'TickDir','out','Layer','top')
xlabel(ax1,'Lon (^oE)','FontSize',13)
ylabel(ax1,'Lat (^oN)','FontSize',13)

title(ax1, [ 'SWOT DA: ' datestr(t_swotda,'dd-mmm HH:MM')], ...
    'FontWeight','bold','FontSize',13)
text(0.82,40.65,'(a)','FontWeight','bold','FontSize',13)

%% ------------------------------------------------------------
% (b) Joint DA
%% ------------------------------------------------------------
ax2 = nexttile(tl,2); hold(ax2,'on')

hp2 = pcolor(ax2, data_joint.lon, data_joint.lat, zeta_joint_anom*100);
shading(ax2,'flat')

[c2,b2] =contour(ax2, data_joint.lon, data_joint.lat, zeta_joint_anom*100, ...
    cont_lev, 'k', 'LineWidth', 0.8);
clabel(c2,b2,'fontsize',14)

plot(ax2, ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)
plot(ax2, poly_lon, poly_lat, '-', 'Color', [1 1 1], 'LineWidth', 2.0)

% ===== SWOT TRACKS USED FOR MASK =====
L1_lon = data_swot.longitude(39+1,:);
L1_lat = data_swot.latitude(39+1,:);
L2_lon = data_swot.longitude(end-4,:);
L2_lat = data_swot.latitude(end-4,:);

% ===== BUILD CLOSED POLYGON =====
poly_lon = [L1_lon fliplr(L2_lon) L1_lon(1)];
poly_lat = [L1_lat fliplr(L2_lat) L1_lat(1)];

% ===== PLOT CLOSED POLYGON =====
plot(ax2,poly_lon, poly_lat, '-', 'Color', colors(4,:), 'LineWidth', 2.2)

% ===== CREATE MASK =====
[Lon2D, Lat2D] = meshgrid(data_joint.lon, data_joint.lat);
in = inpolygon(Lon2D, Lat2D, poly_lon, poly_lat);

alpha_mask = 0.35 * ones(size(in));
alpha_mask(in) = 1;

% ===== APPLY TRANSPARENCY =====
set(hp2,'AlphaData',alpha_mask);
set(hp2,'FaceAlpha','flat');
set(hp2,'AlphaDataMapping','none');

% ADCP transect 2: solid
% plot(ax2, adcp2_lon, adcp2_lat, 'w-', 'LineWidth', 3.0)
plot(ax2, adcp2_lon, adcp2_lat, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2)

% CTD stations
plot(ax2, CTD_lon, CTD_lat,'o','MarkerFaceColor',[0.5210    0.0860    0.8190],'MarkerSize',6)
% Glider track
plot(ax2, glider1.longitude(g1_ind), glider1.latitude(g1_ind),'b-', 'LineWidth', 2)
plot(ax2, glider2.longitude(g2_ind), glider2.latitude(g2_ind),'b-', 'LineWidth', 2)

axis(ax2,'equal')
axis(ax2, axis_def)
box(ax2,'on')
caxis(ax2, cax_lim)

set(ax2,'FontSize',12,'TickDir','out','Layer','top')
xlabel(ax2,'Lon (^oE)','FontSize',13)
set(ax2,'YTickLabel',[])

title(ax2, [ 'Joint DA : ' datestr(t_swotda,'dd-mmm HH:MM')], ...
    'FontWeight','bold','FontSize',13)
text(0.82,40.65,'(b)','FontWeight','bold','FontSize',13)

%% ------------------------------------------------------------
% (c) SWOT
%% ------------------------------------------------------------
ax3 = nexttile(tl,3); hold(ax3,'on')

pcolor(ax3, data_swot.longitude(36:end,:), data_swot.latitude(36:end,:), ...
    swot_ssh(36:end,:));
shading(ax3,'flat')

[c3,b3] = contour(ax3, data_swot.longitude(36:end,:), data_swot.latitude(36:end,:), ...
    swot_ssh(36:end,:), cont_lev, 'k', 'LineWidth', 0.8);
clabel(c3,b3,'fontsize',14)

plot(ax3, ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

% ADCP transect 2: solid
% plot(ax3, adcp2_lon, adcp2_lat, 'w-', 'LineWidth', 3.0)
p1 = plot(ax3, adcp2_lon, adcp2_lat, '-', 'Color', [0.85 0.10 0.10], 'LineWidth',2);
% CTD stations
p2 = plot(ax3, CTD_lon, CTD_lat,'o','MarkerFaceColor',[0.5210    0.0860    0.8190],'MarkerSize',6);
p3 = plot(ax3, glider1.longitude(g1_ind), glider1.latitude(g1_ind),'b-', 'LineWidth', 2);
plot(ax3, glider2.longitude(g2_ind), glider2.latitude(g2_ind),'b-', 'LineWidth', 2)


axis(ax3,'equal')
axis(ax3, axis_def)
box(ax3,'on')
caxis(ax3, cax_lim)

set(ax3,'FontSize',12,'TickDir','out','Layer','top')
xlabel(ax3,'Lon (^oE)','FontSize',13)
set(ax3,'YTickLabel',[])

title(ax3, [ 'SWOT: ' datestr(t_swotda,'dd-mmm HH:MM')], ...
    'FontWeight','bold','FontSize',13)
text(0.82,40.65,'(c)','FontWeight','bold','FontSize',13)

% Reference latitude (pick central latitude of your domain)
lat_ref = 39.5;   % around your figure's center
km = 20;          % length of scale bar

% Convert 20 km into degrees
deg_lon = km / (111 * cosd(lat_ref)); % longitude degrees
deg_lat = km / 111;                   % latitude degrees
lon0 = 2.65;
lat0 = 39.65;              
               
                
hold on
% Draw east-west bar of 20 km
plot([0.86 0.86+deg_lon],[39.7 39.7],'k-','LineWidth',3)
% Add text label above bar
text(0.97,39.76,'20 km','HorizontalAlignment','center','FontSize',14)
text(lon0+deg_lon/2,lat0,'Mallorca','HorizontalAlignment','center','FontSize',16)

legend([p1,p3,p2],{'ADCP','Glider','CTD'},'Box','off','Location','southwest')


%% ============================================================
% Shared colorbar
%% ============================================================
% cb = colorbar(ax3,'eastoutside');
% cb.Label.String = 'SSH anomaly (cm)';
% cb.FontSize = 12;
% cb.TickDirection = 'out';

%% ============================================================
% Optional export
%% ============================================================
print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_SWOT_plots', '-dpng', '-r300');