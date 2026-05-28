% ============================================================
% Figure 2 release script
% Panels: Figure 2a-p
% Purpose: 2D SSH maps on 6 May and assimilation-period SSH evolution/RMSE.
%
% Notes for release:
% - Set plot_index at the top of the script as in the original workflow.
% - plot_index = 5 generates the revised spatial SSH comparison used for the main figure.
% - Other plot_index options support the temporal-evolution and RMSE diagnostics used to assemble Figure 2.
% - Keep figure settings unchanged unless you are intentionally updating the published figure.
% ============================================================

clearvars; clc; % close all
addpath cmocean
addpath TEOS_10
addpath TEOS_10/library/
addpath export_fig-master
RT = load("rt_colormaps.mat");
load('python_colormap_spectral_r.mat');

%% Plot option selector

%   2: temporal evolution of Joint DA performance
%   3: full assimilation-period RMSE time series
%   5: spatial SSH comparison across simulations
% Run only one option at a time to limit memory use.
% Select one option below.
plot_index = 3; % options: 2, 3, or 5

%% DA data

WMOP_path_v1 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v6/'; % GEN
WMOP_path_v3 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v7/'; % SWOT DA
WMOP_path_v4 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v8/'; % In-situ DA
WMOP_path_v5 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9/'; % Joint DA

%% Model file lists

load('coastline_full_westmed_nolakes.mat')

% Get list of files for each version
dir_WMOP_v1 = dir(fullfile([WMOP_path_v1 '*.nc']));
dir_WMOP_v3 = dir(fullfile([WMOP_path_v3 '*.nc']));
dir_WMOP_v4 = dir(fullfile([WMOP_path_v4 '*.nc']));
dir_WMOP_v5 = dir(fullfile([WMOP_path_v5 '*.nc']));

version_labels = { ... % model experiment labels
    'GEN', ...
    'SWOT DA', ...
    'In-situ DA', ...
    'Joint DA'};

% [~, data_v1] = read_nc_file_struct([WMOP_path_v1 filename_WMOP]);
colors =[ 0.0660    0.4430    0.7450
          0.2310    0.6660    0.1960
          0.5210    0.0860    0.8190
          0.9608    0.4667    0.1608
          ];

    axis_def      =  [0.8 2.2 39.2 40.7];
    axis_def_mask =  [0.8 3 39 41];
    SWOT_path     = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';

%% Main processing
count = 1;
% Create a struct array to hold all versions
dir_WMOP_all = {dir_WMOP_v1, dir_WMOP_v3, dir_WMOP_v4, dir_WMOP_v5}; %, dir_WMOP_v6};
WMOP_paths   = {WMOP_path_v1, WMOP_path_v3, WMOP_path_v4, WMOP_path_v5}; %, WMOP_path_v6};

%% ADD ADCP and glider
[~, glider1] = read_nc_file_struct(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/dep0007_sdeep09_scb-sldeep009_L2_2023-04-25_data_dt.nc']);
glider1.Time = datenum(seconds(glider1.time) + datenum(1970,01,01));
% ind_g1 = find(glider1.Time > datenum(2023,5,2,15,43,00) & glider1.Time < datenum(2023,5,4,15,40,00));
% ind_g1 = find(glider1.Time > datenum(2023,4,28,07,58,00) & glider1.Time < datenum(2023,5,2,15,43,00));
ind_g1 = find(glider1.Time > datenum(2023,5,3,15,39,00) & glider1.Time < datenum(2023,5,7,07,38,00));

adcp_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/ADCP/';
[~, data_Leg] = read_nc_file_struct([adcp_path 'ADCP_OS150_FaSt-SWOT_LEG1.nc']);
data_Leg.Time = data_Leg.time + datenum(2023,1,1);
% Select the validation transect.
ind_adcp = find(data_Leg.Time > datenum(2023,04,26,12,00,00) &...
data_Leg.Time < datenum(2023,04,26,18,00,00));

fig4 = figure('Units','pixels','Position',[100 100 1600 900]);

nrow = 2;
ncol = 5;

gap_x = 0.02;
gap_y = 0.06;

margin_left  = 0.06;
margin_right = 0.04;
margin_top   = 0.06;
margin_bot   = 0.07;

panel_width  = (1-margin_left-margin_right-(ncol-1)*gap_x)/ncol;
panel_height = (1-margin_top-margin_bot-(nrow-1)*gap_y)/nrow;

ax = gobjects(nrow*ncol,1);

if plot_index == 5

    for ii = 16   % For making Fig 2 %1:22 %12: 2 May; 19: 9 May
    %% for ii = 1 : 16 % For getting RMSE for table 1

        % === Read file info ===
        filename_WMOP = dir_WMOP_v1(ii).name;
        ind_name = filename_WMOP(18:25);

        [~, data_v1] = read_nc_file_struct([WMOP_path_v1 filename_WMOP]);
        data_v1.time = datenum(seconds(data_v1.ocean_time) + datenum(1968,5,23));
        data_WMOP_all{1} = data_v1;

        dir_SWOT = dir(fullfile([SWOT_path '/*.nc']));
        SWOT_path_new = [SWOT_path '/'];
        ind_name_swot = filename_WMOP(18:25);

        for k = 1:length(dir_SWOT) % for model
            if contains(dir_SWOT(k).name, ind_name_swot)
                [~, data_swot] = read_nc_file_struct([SWOT_path_new dir_SWOT(k).name]);

                data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));

                % Target point
                lon_target = 2;
                lat_target = 40;

                % Compute distance to target point at each grid point and time step
                dist = sqrt((data_swot.longitude - lon_target).^2 + (data_swot.latitude - lat_target).^2);

                % For each time step, find the closest grid point
                [min_dist_per_time, ~] = min(dist, [], 1);

                % Find the time index with the smallest distance overall
                [~, best_time_idx] = min(min_dist_per_time);

                % Get the corresponding time
                data_swot.closest_time = data_swot.TIME(best_time_idx);

                found = true;
                break;
            end
        end

        % Read v3-v6
        for ver = 1:4 %2:5
            found = false;
            this_dir = dir_WMOP_all{ver};
            this_path = WMOP_paths{ver};

            for k = 1:length(this_dir) % for model
                if contains(this_dir(k).name, ind_name)
                    [~, data_tmp] = read_nc_file_struct([this_path this_dir(k).name]);
                    data_tmp.time = datenum(seconds(data_tmp.ocean_time) + datenum(1968,5,23));
                    data_WMOP_all{ver} = data_tmp;
                    found = true;
                    break;
                end
            end

            if ~found
                fprintf('No file found for %s in v%d\n', ind_name, ver + (ver > 1));
                data_WMOP_all{ver} = [];
            end
        end

        % === Frame Loop ===
        Nt = size(data_v1.zeta, 3);

        %% Hour-selection block
        j_val = 6; % 7;

        for j = j_val %1: 2 : Nt

            clf(fig4);
            set(fig4,'Color','w','Position',[80 120 1350 950]);
            tl = tiledlayout(fig4,2,3,'TileSpacing','compact','Padding','compact');

            panel_letters = {'(a)','(b)','(c)','(d)'};

fig_loc = [1,2,4,5];
            for k = 1:4

                nexttile(fig_loc(k))

                clear zeta
                data_WMOP = data_WMOP_all{k};

                if ~isempty(data_WMOP) && isfield(data_WMOP, 'zeta') && j <= size(data_WMOP.zeta,3)

                    zeta = data_WMOP.zeta(:,:,j)';

                    lon_min = axis_def_mask(1);
                    lon_max = axis_def_mask(2);
                    lat_min = axis_def_mask(3);
                    lat_max = axis_def_mask(4);
                    [Lon, Lat] = meshgrid(data_WMOP.lon, data_WMOP.lat);

                    % Create mask for the desired region
                    mask = Lon >= lon_min & Lon <= lon_max  & ...
                           Lat  >= lat_min & Lat  <= lat_max ;

                    zeta_masked = zeta;
                    zeta_masked(~mask) = NaN;

                    zeta_anom = zeta - nanmean(zeta_masked(:));

                    % ===== SSH FIELD =====
                    hp = pcolor(data_WMOP.lon, data_WMOP.lat, zeta_anom*100);
                    shading flat
                    hold on

                    % ===== CONTOURS =====
                    [C1,h1] = contour(data_WMOP.lon, data_WMOP.lat, zeta_anom*100, -10:1:10, ...
                            'k','LineWidth',1.2);
                    clabel(C1,h1,'FontSize',12,'Color','k')
                    hold on

                    % ===== SWOT TRACKS USED FOR MASK =====
                    L1_lon = data_swot.longitude(39+1,:);
                    L1_lat = data_swot.latitude(39+1,:);
                    L2_lon = data_swot.longitude(end-4,:);
                    L2_lat = data_swot.latitude(end-4,:);

                    % ===== BUILD CLOSED POLYGON =====
                    poly_lon = [L1_lon fliplr(L2_lon) L1_lon(1)];
                    poly_lat = [L1_lat fliplr(L2_lat) L1_lat(1)];

                    % ===== PLOT CLOSED POLYGON =====
                    plot(poly_lon, poly_lat, '-', 'Color', colors(k,:), 'LineWidth', 2.2)

                    % ===== CREATE MASK =====
                    [Lon2D, Lat2D] = meshgrid(data_WMOP.lon, data_WMOP.lat);
                    in = inpolygon(Lon2D, Lat2D, poly_lon, poly_lat);

                    alpha_mask = 0.35 * ones(size(in));
                    alpha_mask(in) = 1;

                    % ===== APPLY TRANSPARENCY =====
                    set(hp,'AlphaData',alpha_mask);
                    set(hp,'FaceAlpha','flat');
                    set(hp,'AlphaDataMapping','none');

                    % ===== COASTLINE =====
                    plot(ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

                    % ===== AXIS SETTINGS =====
                    axis equal
                    axis(axis_def)
                    box on
                    colormap(gca, RT.rt_colormaps.section)
                    caxis([-7 7])

                    if k == 1
                        set(gca,'FontSize',14,'TickDir','out','Layer','top','YTick',39:0.5:41,'XTickLabel',[])
                        ylabel('Lat (^oN)')
                    elseif k == 2
                        set(gca,'FontSize',14,'TickDir','out','Layer','top','YTick',39:0.5:41,'XTickLabel',[],'YTickLabel',[])
                    elseif k == 3
                        set(gca,'FontSize',14,'TickDir','out','Layer','top','YTick',39:0.5:41)
                        ylabel('Lat (^oN)')
                        xlabel('Lon (^oE)')
                    elseif k == 4
                        set(gca,'FontSize',14,'TickDir','out','Layer','top','YTick',39:0.5:41,'YTickLabel',[])
                        xlabel('Lon (^oE)')
                    end

                    if k == 3
                        xlabel('Lon (^oE)')
                    % else
                    %     set(gca,'XTickLabel',[])
                    end

                    title([panel_letters{k} ' ' version_labels{k} ' 06 May'], ...
                          'FontWeight','bold','FontSize',13)

                    %% ---------
                    % transect_line = [1.33485 40.3721; 1.74528 39.1444];
                    transect_line = [1.48721 40.1672; 1.68378 39.5468];

                    lon_target = linspace(transect_line(1,1), transect_line(2,1), 100);
                    lat_target = linspace(transect_line(1,2), transect_line(2,2), 100);
                    ind_S = find(lat_target < 39.2);
                    lon_target(ind_S) = NaN;
                    lat_target(ind_S) = NaN;

                    hold on
                    p1 = plot(glider1.longitude(ind_g1), glider1.latitude(ind_g1), 'b-', 'LineWidth', 2);

                    % Reference latitude for distance conversion
                    lat_ref = 39.5;
                    km = 20;

                    % Convert 20 km into degrees
                    deg_lon = km / (111 * cosd(lat_ref));
                    deg_lat = km / 111;

                    % Choose where to place the bar (lower left corner)
                    lon0 = 2.65;
                    lat0 = 39.65;

                end

                if k == 1
                    %% SWOT
                    lon_min = axis_def_mask(1);
                    lon_max = axis_def_mask(2);
                    lat_min = axis_def_mask(3);
                    lat_max = axis_def_mask(4);

                    mask = data_swot.longitude >= lon_min  & data_swot.longitude <= lon_max  & ...
                           data_swot.latitude  >= lat_min  & data_swot.latitude  <= lat_max ;

                    ssha_masked   = data_swot.ssha_filtered;
                    swot_lon_mask = data_swot.longitude;
                    swot_lat_mask = data_swot.latitude;

                    ssha_masked(~mask)   = NaN;
                    swot_lon_mask(~mask) = NaN;
                    swot_lat_mask(~mask) = NaN;

                    swot_ssh = (ssha_masked - nanmean(nanmean(ssha_masked))) * 100;
                end

                model_to_swot = interp2(data_WMOP.lon, data_WMOP.lat, zeta_anom*100, ...
                                        swot_lon_mask(36:end,:), swot_lat_mask(36:end,:));
                A = swot_ssh(36:end,:);
                B = model_to_swot;
                RMSE = sqrt(nanmean((A(:) - B(:)).^2));

                RMSE_all(ii,k) = RMSE;

                %% --- RMSE label in upper-right corner ---
                text(0.97, 0.95, ['RMSE ' num2str(round(RMSE,2)) ' cm'], ...
                    'Units','normalized', ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','top', ...
                    'FontSize',12, ...
                    'FontWeight','bold', ...
                    'Color',colors(k,:), ...
                    'BackgroundColor','w', ...
                    'Margin',3);

            end

            %% ===============================
            % PANEL (e): SWOT
            %% ===============================
            nexttile(3)

            swot_transect_line = [1.48721 40.1672; 1.68378 39.5468];
            ssha = (ssha_masked - nanmean(nanmean(ssha_masked))) * 100;

            pcolor(data_swot.longitude(36:end,:), data_swot.latitude(36:end,:), ssha(36:end,:));
            shading interp
            hold on

            axis equal
            axis(axis_def)

            % === Add Coastline ===
            plot(ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

            hold on
            [C2,h2] = contour(data_swot.longitude(36:end,:), data_swot.latitude(36:end,:), ...
                    ssha(36:end,:), -10:1:10, 'k', 'LineWidth', 1.2)
            clabel(C2,h2,'FontSize',12,'Color','k')
            caxis([-7 7])
            colormap(gca, RT.rt_colormaps.section)

            set(gca, 'FontSize',14, 'Layer','top', 'TickDir','out', ...
                     'YTick',39:0.5:41, 'YTickLabel', []);
            box on

            p1 = plot(glider1.longitude(ind_g1), glider1.latitude(ind_g1), 'b-', 'LineWidth', 2);

            title('(e) SWOT', 'FontWeight','bold', 'FontSize',13)
            xlabel('Lon (^oE)')

% Reference latitude for distance conversion
lat_ref = 39.5;   % central latitude for the analysis domain
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

        end

        print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig2_upper_v2'], '-dpng', '-r300');

    end
end

%% Joint DA temporal-evolution settings
ti = [8 8 7 7 6];
if plot_index == 2 % Joint DA temporal evolution and collocated SWOT observations

%% Joint DA and SWOT comparison panels
        %% ===============================
        % ROW 2: JOINT DA EVOLUTION
        %% ===============================
        joint_dates = [ ...
            datenum(2023,4,23)
            datenum(2023,4,25)
            datenum(2023,4,27)
            datenum(2023,4,30)
            datenum(2023,5,4)];

    titles_top = { ...
        '(f) Joint DA 23 Apr', ...
        '(g) Joint DA 25 Apr', ...
        '(h) Joint DA 27 Apr', ...
        '(i) Joint DA 30 Apr', ...
        '(j) Joint DA 04 May'};

    titles_bot = { ...
        '(k) SWOT 23 Apr', ...
        '(l) SWOT 25 Apr', ...
        '(m) SWOT 27 Apr', ...
        '(n) SWOT 30 Apr', ...
        '(o) SWOT 04 May'};

% manually assign Joint DA files for each date
mod = 0;
joint_file_ids = [3 5 7 10 14]; %[1+mod 4+mod 7+mod 10+mod 13+mod];
SWOT_ind       = [21 23 25 28 32]; %[19+mod, 22+mod, 25+mod, 28+mod, 31+mod];

figure('Color','w','Position',[80 80 1700 760]) %,'Visible','off');
tl = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');

        for i = 1:5

    %% --- load Joint DA file assigned to this date ---
    file_id = joint_file_ids(i);
    [~, data_joint] = read_nc_file_struct([WMOP_path_v5 dir_WMOP_v5(file_id).name]);
    data_joint.time = datenum(seconds(data_joint.ocean_time) + datenum(1968,5,23));

     % [~,ti] = min(abs(data_joint.time - joint_dates(i)));

    zeta_bot = data_joint.zeta(:,:,ti(i))'; % selected hourly snapshot

            [Lon_bot, Lat_bot] = meshgrid(data_joint.lon, data_joint.lat);
            mask_bot = Lon_bot >= axis_def_mask(1) & Lon_bot <= axis_def_mask(2) & ...
                       Lat_bot >= axis_def_mask(3) & Lat_bot <= axis_def_mask(4);

            zeta_bot_masked = zeta_bot;
            zeta_bot_masked(~mask_bot) = NaN;

            zeta_anom_bot = zeta_bot - nanmean(zeta_bot_masked(:));

            SWOT_path_new = [SWOT_path '/'];
            dir_SWOT    = dir(fullfile([SWOT_path  '/*.nc']));

        ind_name =  dir_WMOP_v5(file_id).name(18:25);

        for k = 1:length(dir_SWOT) % for model
            if contains(dir_SWOT(k).name, ind_name)
                [info, data_swot_bot] = read_nc_file_struct([SWOT_path_new dir_SWOT(k).name]);
            end
        end

            % [~, data_swot_bot] = read_nc_file_struct([SWOT_path_new dir_SWOT(SWOT_ind(i)).name]);
            data_swot_bot.TIME = datenum(seconds(data_swot_bot.time) + datenum(2000,1,1));

        lon_target = 2;
        lat_target = 40;

        dist = sqrt((data_swot_bot.longitude - lon_target).^2 + ...
                    (data_swot_bot.latitude  - lat_target).^2);

        [min_dist_per_time, ~] = min(dist, [], 1);
        [~, best_time_idx] = min(min_dist_per_time);

        swot_time_used(i) = data_swot_bot.TIME(best_time_idx);

    %% --- SWOT mask in selected domain ---
    lon_min = axis_def_mask(1);
    lon_max = axis_def_mask(2);
    lat_min = axis_def_mask(3);
    lat_max = axis_def_mask(4);

    mask_swot_bot = data_swot_bot.longitude >= lon_min  & data_swot_bot.longitude <= lon_max  & ...
                    data_swot_bot.latitude  >= lat_min  & data_swot_bot.latitude  <= lat_max;

    ssha_masked_bot   = data_swot_bot.ssha_filtered;
    swot_lon_mask_bot = data_swot_bot.longitude;
    swot_lat_mask_bot = data_swot_bot.latitude;

    ssha_masked_bot(~mask_swot_bot)   = NaN;
    swot_lon_mask_bot(~mask_swot_bot) = NaN;
    swot_lat_mask_bot(~mask_swot_bot) = NaN;

    swot_ssh_bot = (ssha_masked_bot - nanmean(ssha_masked_bot(:)))*100;

            % add same SWOT mask polygon for consistency
            L1_lon = data_swot_bot.longitude(39+1,:);
            L1_lat = data_swot_bot.latitude(39+1,:);
            L2_lon = data_swot_bot.longitude(end-4,:);
            L2_lat = data_swot_bot.latitude(end-4,:);

            poly_lon = [L1_lon fliplr(L2_lon) L1_lon(1)];
            poly_lat = [L1_lat fliplr(L2_lat) L1_lat(1)];

            [Lon2D_bot, Lat2D_bot] = meshgrid(data_joint.lon, data_joint.lat);
            in_bot = inpolygon(Lon2D_bot, Lat2D_bot, poly_lon, poly_lat);

            alpha_mask_bot = 0.45 * ones(size(in_bot));
            alpha_mask_bot(in_bot) = 1;

        %% --- RMSE between model and SWOT over SWOT swath ---
        model_to_swot_bot = interp2(data_joint.lon, data_joint.lat, zeta_anom_bot*100, ...
                                    swot_lon_mask_bot(36:end,:), swot_lat_mask_bot(36:end,:));

        A_bot = swot_ssh_bot(36:end,:);
        B_bot = model_to_swot_bot;

        RMSE_bot(i) = sqrt(nanmean((A_bot(:) - B_bot(:)).^2));

 %% ===============================
        % UPPER ROW: JOINT DA (f-j)
        %% ===============================
        ax(i) = nexttile(i);

        hp2 = pcolor(data_joint.lon, data_joint.lat, zeta_anom_bot*100);
        shading flat
        hold on

        [C2,h2] = contour(data_joint.lon, data_joint.lat, zeta_anom_bot*100, -10:1:10, ...
                          'k','LineWidth',1.3);
        clabel(C2,h2,'FontSize',12,'Color','k')

        plot(ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

        % show SWOT swath outline only, no SWOT contours
        plot(poly_lon, poly_lat, '-', 'Color', colors(4,:), 'LineWidth', 2.2)

        set(hp2,'AlphaData',alpha_mask_bot);
        set(hp2,'FaceAlpha','flat');
        set(hp2,'AlphaDataMapping','none');

        axis equal
        axis(axis_def)
        caxis([-7 7])
        colormap(ax(i), RT.rt_colormaps.section)

        if i == 1
            ylabel('Lat (^oN)')
            set(gca,'YTick',39:0.5:41,'FontSize',14,'Layer','top','TickDir','out')
        else
            set(gca,'YTick',39:0.5:41,'YTickLabel',[],'FontSize',14,'Layer','top','TickDir','out')
        end

        set(gca,'XTickLabel',[])
        box on
        title(titles_top{i}, 'FontWeight','bold', 'FontSize',14)

        text(0.97, 0.95, ['RMSE ' num2str(round(RMSE_bot(i),2)) ' cm'], ...
            'Units','normalized', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','top', ...
            'FontSize',12, ...
            'FontWeight','bold', ...
            'Color',colors(4,:), ...
            'BackgroundColor','w', ...
            'Margin',3);

        %% ===============================
        % LOWER ROW: SWOT (k-o)
        %% ===============================
        ax(i+5) = nexttile(i+5);

        hp3 = pcolor(data_swot_bot.longitude(36:end,:), data_swot_bot.latitude(36:end,:), swot_ssh_bot(36:end,:));
        shading flat
        hold on

        [Cs, hs] = contour(data_swot_bot.longitude(36:end,:), data_swot_bot.latitude(36:end,:), swot_ssh_bot(36:end,:), ...
                           -10:1:10, 'k', 'LineWidth', 1.3);
        clabel(Cs, hs, 'FontSize',12, 'Color','k')

        plot(ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.0)

        axis equal
        axis(axis_def)
        caxis([-7 7])
        colormap(ax(i+5), RT.rt_colormaps.section)

        if i == 1
            ylabel('Lat (^oN)')
            set(gca,'YTick',39:0.5:41,'FontSize',14,'Layer','top','TickDir','out')
        else
            set(gca,'YTick',39:0.5:41,'YTickLabel',[],'FontSize',14,'Layer','top','TickDir','out')
        end

        if i == 3
            xlabel('Lon (^oE)')
        end

        box on
        title(titles_bot{i}, 'FontWeight','bold', 'FontSize',14)

        end
        %% End Joint DA and SWOT comparison panels
        % print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig2_middle'], '-dpng', '-r300');

end

% [info, AAA] = read_nc_file_struct([SWOT_path dir_SWOT(20).name]);

if plot_index == 3

    joint_file_ids = 1:16;
    RMSE_bot = nan(length(joint_file_ids),3);
    swot_time_used = nan(length(joint_file_ids),1);
    model_time_used = nan(length(joint_file_ids),3);

    dir_SWOT = dir(fullfile([SWOT_path '*.nc']));

    for i = 1:length(joint_file_ids)
        clear data_swot_bot data_joint %RMSE_bot
        file_id = joint_file_ids(i);

        %% ----------------------------------------------------
        % First find corresponding SWOT file and its actual time
        %% ----------------------------------------------------
        ind_name = dir_WMOP_v5(file_id).name(18:25);

        data_swot_bot = [];
        for k = 1:length(dir_SWOT)
            if contains(dir_SWOT(k).name, ind_name)
                [~, data_swot_bot] = read_nc_file_struct([SWOT_path dir_SWOT(k).name]);
                break
            end
        end

        if isempty(data_swot_bot)
            warning('No SWOT file found for %s', ind_name)
            continue
        end

        data_swot_bot.TIME = datenum(seconds(data_swot_bot.time) + datenum(2000,1,1));

        lon_target = 2;
        lat_target = 40;

        dist = sqrt((data_swot_bot.longitude - lon_target).^2 + ...
                    (data_swot_bot.latitude  - lat_target).^2);

        [min_dist_per_time, ~] = min(dist, [], 1);
        [~, best_time_idx] = min(min_dist_per_time);

        swot_time_used(i) = data_swot_bot.TIME(best_time_idx);

        %% ----------------------------------------------------
        % SWOT mask and anomaly
        %% ----------------------------------------------------
        lon_min = axis_def_mask(1);
        lon_max = axis_def_mask(2);
        lat_min = axis_def_mask(3);
        lat_max = axis_def_mask(4);

        mask_swot_bot = data_swot_bot.longitude >= lon_min & data_swot_bot.longitude <= lon_max & ...
                        data_swot_bot.latitude  >= lat_min & data_swot_bot.latitude  <= lat_max;

        ssha_masked_bot   = data_swot_bot.ssha_filtered;
        swot_lon_mask_bot = data_swot_bot.longitude;
        swot_lat_mask_bot = data_swot_bot.latitude;

        ssha_masked_bot(~mask_swot_bot)   = NaN;
        swot_lon_mask_bot(~mask_swot_bot) = NaN;
        swot_lat_mask_bot(~mask_swot_bot) = NaN;

        swot_ssh_bot = (ssha_masked_bot - nanmean(ssha_masked_bot(:))) * 100;

        %% ----------------------------------------------------
        % Loop over the 3 DA experiments
        %% ----------------------------------------------------
        for j = 1:4
            if j == 1
                [~, data_joint] = read_nc_file_struct([WMOP_path_v1 dir_WMOP_v1(file_id).name]); % v7
            elseif j == 2
                [~, data_joint] = read_nc_file_struct([WMOP_path_v3 dir_WMOP_v3(file_id).name]); % v7
            elseif j == 3
                [~, data_joint] = read_nc_file_struct([WMOP_path_v4 dir_WMOP_v4(file_id).name]); % v8
            elseif j == 4
                [~, data_joint] = read_nc_file_struct([WMOP_path_v5 dir_WMOP_v5(file_id).name]); % v9
            end

            data_joint.time = datenum(seconds(data_joint.ocean_time) + datenum(1968,5,23));

            % use closest model time to SWOT observation time
            [~, ti] = min(abs(data_joint.time - swot_time_used(i)));
            model_time_used(i,j) = data_joint.time(ti);

            zeta_bot = data_joint.zeta(:,:,ti)';

            [Lon_bot, Lat_bot] = meshgrid(data_joint.lon, data_joint.lat);
            mask_bot = Lon_bot >= axis_def_mask(1) & Lon_bot <= axis_def_mask(2) & ...
                       Lat_bot >= axis_def_mask(3) & Lat_bot <= axis_def_mask(4);

            zeta_bot_masked = zeta_bot;
            zeta_bot_masked(~mask_bot) = NaN;

            zeta_anom_bot = zeta_bot - nanmean(zeta_bot_masked(:));

            %% --- RMSE between model and SWOT over SWOT swath ---
            model_to_swot_bot = interp2(data_joint.lon, data_joint.lat, zeta_anom_bot*100, ...
                                        swot_lon_mask_bot(36:end,:), swot_lat_mask_bot(36:end,:));

            A_bot = swot_ssh_bot(36:end,:);
            B_bot = model_to_swot_bot;

            RMSE_bot(i,j) = sqrt(nanmean((A_bot(:) - B_bot(:)).^2));
        end
    end

    %% --------------------------------------------------------
    % Sort by actual SWOT time
    %% --------------------------------------------------------
    [swot_time_used, isort] = sort(swot_time_used);
    RMSE_bot = RMSE_bot(isort,:);
    model_time_used = model_time_used(isort,:);

    tplot = datetime(swot_time_used,'ConvertFrom','datenum');

    %% --------------------------------------------------------
    % Plot

    %% --------------------------------------------------------
    figure('color','w','position',[200 200 1050 340]); hold on

    % set(fig4,'Color','w','Position',[80 120 1350 950]);

    p0 = plot(tplot, RMSE_bot(:,1), 'o-', 'LineWidth',2.2, ...
        'Color',[0.0660    0.4430    0.7450], 'MarkerFaceColor',[0.0660    0.4430    0.7450],'MarkerSize', 7);

    p1 = plot(tplot, RMSE_bot(:,2), 'o-', 'LineWidth',2.2, ...
        'Color',[0.2310 0.6660 0.1960],'MarkerFaceColor',[0.2310 0.6660 0.1960],'MarkerSize', 7);

    p2 = plot(tplot, RMSE_bot(:,3), 'o-', 'LineWidth',2.2, ...
        'Color',[0.5210 0.0860 0.8190],'MarkerFaceColor',[0.5210 0.0860 0.8190],'MarkerSize', 7);

    p3 = plot(tplot, RMSE_bot(:,4), 'o-', 'LineWidth',2.2, ...
        'Color',[0.9608 0.4667 0.1608],'MarkerFaceColor',[0.9608 0.4667 0.1608],'MarkerSize', 7);
%1:3:18
    % plot(tplot([3,5,7,10,14,16]), RMSE_bot([3,5,7,10,14,16],3), 'o', ...
    %     'Color',[0.9608 0.4667 0.1608], ...
    %     'MarkerSize',8, ...
    %     'MarkerFaceColor',[0.9608 0.4667 0.1608], ...
    %     'MarkerEdgeColor','k')

    % plot([datetime(2023,4,26,00,00,00) datetime(2023,4,26,00,00,00)], [1.5 4.5], 'r--')

    box on
    set(gca,'FontSize',14,'LineWidth',1.5,...
            'TickDir','out',...
            'FontName','Helvetica')

    ylabel('RMSE (cm)')
    xlabel('Date')

    xtickformat('dd-MMM')
    % xtickangle(30)
    xlim([datetime(2023,4,21,00,00,00) tplot(end)])

    legend([p0,p1,p2,p3], {'GEN','SWOT DA','In-situ DA','Joint DA'}, ...
           'Location','northeast', ...
           'Box','off')

    text(tplot(1), 4.3, '(f)', 'FontWeight','bold', 'FontSize',14)

    grid on
    ax = gca;
    ax.GridAlpha = 0.15;
    ylim([1.4 4.5])

    print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig2_lower_v2'], '-dpng', '-r300');

    % disp(table(tplot(:), RMSE_bot(:,1), RMSE_bot(:,2), RMSE_bot(:,3), ...
    %     datetime(model_time_used(:,1),'ConvertFrom','datenum'), ...
    %     datetime(model_time_used(:,2),'ConvertFrom','datenum'), ...
    %     datetime(model_time_used(:,3),'ConvertFrom','datenum'), ...
    %     'VariableNames', {'SWOT_Time','SWOT_DA','Insitu_DA','Joint_DA', ...
    %                       'ModelTime_v7','ModelTime_v8','ModelTime_v9'}))
end
