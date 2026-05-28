% ============================================================
% Figure 1 release script
% Panels: Figure 1a-d
% Purpose: Observational context and data assimilation framework.
%
% Notes for release:
% - Run this script directly.
% - Requires CTD, ADCP, glider, SWOT, coastline, and colormap assets.
% - Uses local absolute paths as in the submitted analysis version.
% - Keep figure settings unchanged unless you are intentionally updating the published figure.
% ============================================================

clearvars; clc; close all
addpath cmocean
addpath TEOS_10
addpath TEOS_10/library/
addpath Sea-Bird-Toolbox-master/CTD_CNV/
RT = load("rt_colormaps.mat");
load('coastline_full_westmed_nolakes.mat')
load('python_colormap_spectral_r.mat');

axis_def  =  [0.8 3 39 41];

fig = figure('Units','pixels','Position',[100 100 1250 1200]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')


%% CTD
    CTD1_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/CTDs/Leg1/';
    CTD2_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/CTDs/Leg2/';

    CTD_dir = dir(CTD1_path);
    for i = 3 : 14
        cnv = readSBScnv([CTD1_path CTD_dir(i).name]);
        lat1_str = cnv.instrumentheaders.NMEALatitude;   % e.g. '39 41.21 N'
        lon1_str = cnv.instrumentheaders.NMEALongitude;  % e.g. '001 26.46 E'
        tok = regexp(lat1_str,'(\d+)\s+(\d+\.?\d*)\s*([NS])','tokens','once');
        lat1(i,1) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='S', lat1(i,1) = -lat1(i,1); end

        tok = regexp(lon1_str,'(\d+)\s+(\d+\.?\d*)\s*([EW])','tokens','once');
        lon1(i,1) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='W', lon1(i,1) = -lon1(i,1); end
    end

    CTD_dir = dir(CTD2_path);
        for i = 3 : 12
        cnv = readSBScnv([CTD2_path CTD_dir(i).name]);
        lat2_str = cnv.instrumentheaders.NMEALatitude;   % e.g. '39 41.21 N'
        lon2_str = cnv.instrumentheaders.NMEALongitude;  % e.g. '001 26.46 E'
        tok = regexp(lat2_str,'(\d+)\s+(\d+\.?\d*)\s*([NS])','tokens','once');
        lat2(i,1) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='S', lat2(i,1) = -lat2(i,1); end

        tok = regexp(lon2_str,'(\d+)\s+(\d+\.?\d*)\s*([EW])','tokens','once');
        lon2(i,1) = str2double(tok{1}) + str2double(tok{2})/60; if tok{3}=='W', lon2(i,1) = -lon2(i,1); end
    end
    % figure; hold on
    % plot(lon1,lat1,'bo')
    % plot(lon2,lat2,'ko')
    % axis([1 3 38 40])

%% ADCP
    adcp_path = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/ADCP/';

    [~, data_Leg] = read_nc_file_struct([adcp_path 'ADCP_OS150_FaSt-SWOT_LEG2.nc']);
    data_Leg.Time = data_Leg.time + datenum(2023,1,1);
    ind{1} = find(data_Leg.Time > datenum(2023,05,08,02,00,00) &...
        data_Leg.Time < datenum(2023,05,08,9,35,00));

    ind{2} = find(data_Leg.Time > datenum(2023,05,08,11,00,00) &...
        data_Leg.Time < datenum(2023,05,08,19,35,00));

    ind{3} = find(data_Leg.Time > datenum(2023,05,08,20,30,00) &...
        data_Leg.Time < datenum(2023,05,09,03,35,00));

    ind{4} = find(data_Leg.Time > datenum(2023,05,09,04,40,00) &...
        data_Leg.Time < datenum(2023,05,09,14,30,00));

    ind{5} = find(data_Leg.Time > datenum(2023,05,09,15,40,00) &...
        data_Leg.Time < datenum(2023,05,09,23,30,00));

    ind{6} = find(data_Leg.Time > datenum(2023,05,08,02,00,00) &...
        data_Leg.Time < datenum(2023,05,09,23,30,00));

%% Glider
% nexttile(3)
ax3 = axes('Position',[0.08 0.08 0.36 0.36]);

    [~, glider1] = read_nc_file_struct(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/dep0040_sdeep01_scb-sldeep001_L2_2023-04-25_data_dt.nc']);
    glider1.Time = datenum(seconds(glider1.time) + datenum(1970,01,01));

    [info, glider2] = read_nc_file_struct(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/Cruise_data/Glider/dep0007_sdeep09_scb-sldeep009_L2_2023-04-25_data_dt.nc']);
    glider2.Time = datenum(seconds(glider2.time) + datenum(1970,01,01));

    glider2.CT = gsw_CT_from_t(glider2.salinity, glider2.temperature, 0);  % Converts Temperature (t) to Conservative Temperature (θ)
    glider2.SA = gsw_SA_from_SP(glider2.salinity, 0, glider2.longitude, glider2.latitude);  % Converts practical salinity (PSS-78) to Absolute Salinity (SA)
    glider2.PD = gsw_rho(glider2.SA, glider2.CT, 0) -1000;  % In-situ density in kg/m^3
    % glider2.PD = glider2.density -1000;

    glider1.CT = gsw_CT_from_t(glider1.salinity, glider1.temperature, 0);  % Converts Temperature (t) to Conservative Temperature (θ)
    glider1.SA = gsw_SA_from_SP(glider1.salinity, 0, glider1.longitude, glider1.latitude);  % Converts practical salinity (PSS-78) to Absolute Salinity (SA)
    glider1.PD = gsw_rho(glider1.SA, glider1.CT, 0) -1000;  % In-situ density in kg/m^3
    % glider1.PD = glider1.density -1000;

    g1{1} = find(glider1.Time >= datenum(2023,4,28,7,58,00) ...
                    & glider1.Time <= datenum(2023,5,2,15,43,00));

    g1{2} = find(glider1.Time >= datenum(2023,5,2,15,43,00) ...
                    & glider1.Time <= datenum(2023,5,4,15,40,00));

    g1{3} = find(glider1.Time >= datenum(2023,5,4,15,40,00) ...
                    & glider1.Time <= datenum(2023,5,7,23,40,00));

    % g1{4} = find(glider1.Time >= datenum(2023,5,7,23,40,00) ...
    %                 & glider1.Time <= datenum(2023,5,10,07,03,00));

    g2{1} = find(glider2.Time >= datenum(2023,4,28,9,48,00) ...
                    & glider2.Time <= datenum(2023,5,1,16,50,00));

    g2{2} = find(glider2.Time >= datenum(2023,5,1,16,50,00) ...
                    & glider2.Time <= datenum(2023,5,3,15,45,00));

    % g2{3} = find(glider2.Time >= datenum(2023,5,3,15,45,00) ...
    %                 & glider2.Time <= datenum(2023,5,7,07,04,00));

    g2{3} = find(glider2.Time >= datenum(2023,5,3,15,45,00) ...
                & glider2.Time <= datenum(2023,5,6,23,59,59));

    % g1{4} = find(glider2.Time >= datenum(2023,5,7,07,38,00) ...
    %                 & glider2.Time <= datenum(2023,5,10,07,12,00));

    [lat_sorted, idx] = unique(glider2.latitude(g2{3}), 'stable');

    % PDi = glider2.PD(:, g2{3}(idx));
    PDi = glider2.PD(:, g2{3}(idx));
    Si  = glider2.salinity(:, g2{3}(idx));

    pcolor(lat_sorted, -glider2.depth, Si); shading flat
    hold on
    [c, b] = contour(lat_sorted, -glider2.depth, PDi, 28.8:0.05:30.3, 'k','LineWidth',1.5);
    [c1, b1] = contour(lat_sorted, -glider2.depth, PDi, 28.8:0.05:30.3, 'k','LineWidth',1.5);
    clabel(c1, b1, 'Color','k','fontsize',14);

    ylim([-600 -0]);
    % title([ datestr(glider2.Time(g2{3}(1)),  'dd-mmm HH:MM'), ...
    % ' - ', datestr(glider2.Time(g2{3}(end)),'dd-mmm HH:MM') ])
    set(gca,'Layer','top','FontSize',16,'tickdir','out'); box on
    colormap(ax3, RT.rt_colormaps.sst2);
    % cmocean('delta')
    % caxis([28.3 29.1])%
    caxis([38.2 38.6])
    hold on;

    text(39.56,-40,'(b)',...
         'FontWeight','bold','FontSize',16,'Color','w')
    title('Glider Salinity (3 May 17:09 - 6 May 23:38)','Units','normalized',...
         'FontWeight','bold','FontSize',16,'Color','k')
    ylabel('z (m)')
    xlabel('Lat (\circN)');
    % cb = colorbar('position',[0.45 0.08 0.012 0.16]);
    % cb.Label.String = 'g kg^{-1}';
    % cb.YAxisLocation = 'right';
    % set(cb,'FontSize',16)
    % % cb.Yaxis
    % cb.YAxisLocation = 'Left';

% return
%% SWOT
% nexttile(1)
    ax1 = axes('Position',[0.06 0.55 0.40 0.40]);
    % SWOT_path    = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/v2_0/IMEDEA_SWOT_016/';
    SWOT_path    = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/SWOT/L3/SWOT_L3_v3_016/';

    dir_SWOT    = dir(fullfile([SWOT_path '/*.nc']));
    SWOT_path_new = [SWOT_path '/'];

    [~, data_swot] = read_nc_file_struct([SWOT_path_new dir_SWOT(34).name]);
    data_swot.TIME = datenum(seconds(data_swot.time) + datenum(2000,1,1));

    % Unpack axis limits
    lon_min = axis_def(1);
    lon_max = axis_def(2);
    lat_min = axis_def(3);
    lat_max = axis_def(4);
    % Create mask for the desired region
    mask = data_swot.longitude >= lon_min  & data_swot.longitude <= lon_max  & ...
           data_swot.latitude  >= lat_min  & data_swot.latitude  <= lat_max ;

    % Masked data: use NaN outside region
    ssha_masked = data_swot.ssha_filtered;
    ssha_masked(~mask) = NaN;

    pcolor(data_swot.longitude, data_swot.latitude, (ssha_masked - nanmean(nanmean(ssha_masked)))*100); shading interp
    hold on
    axis equal;
    axis(axis_def);
    % === Add Coastline ===
    plot(ncst(:,1), ncst(:,2), 'k', 'LineWidth', 1.2);  % overlay coast in black

    hold on
    [c,b] = contour( data_swot.longitude, data_swot.latitude, (ssha_masked - nanmean(nanmean(ssha_masked)))*100, [-10: 1: 10],'k','LineWidth',2);
    clabel(c,b, 'FontSize',16)
    caxis([-7  7])
    set(gca, 'FontSize',16,'Layer','top','tickdir','out','ytick',[39:0.5:41]); box on

    g7 = find(glider2.Time >= datenum(2023,4,28,9,48,00) ...
                    & glider2.Time <= datenum(2023,5,7,07,04,00));

    % plot(swot_transect_line(:,1), swot_transect_line(:,2),'r-','LineWidth',1.5); %axis equal
    p1 = plot(glider1.longitude(g7), glider1.latitude(g7), 'b-','LineWidth',2);
    % p2 = plot(data_Leg.lon(ind{6}), data_Leg.lat(ind{6}), 'b-','LineWidth',2);

    p2 = plot(lon1,lat1,'o','MarkerFaceColor',[0.5210    0.0860    0.8190],'MarkerSize',8);
    % plot(lon2,lat2,'go','MarkerFaceColor','g')

    % text(0,1.035,'(a) SWOT SSH (3 May 2023 5:34)','Units','normalized','FontWeight','bold','FontSize',16)
    text(0.83,40.85,'(a)','FontWeight','bold','FontSize',16)
    title('SWOT SSH (6 May 2023 5:06)','Units','normalized','FontWeight','bold','FontSize',16)
    xlabel('Lon (\circE)'); ylabel('Lat (\circN)');

    % colors =[ 0.0660    0.4430    0.7450
    %       0.2310    0.6660    0.1960
    %       0.9608    0.4667    0.1608
    %       0.5210    0.0860    0.8190];

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
    plot([1 1+deg_lon],[39.7 39.7],'k-','LineWidth',3)
    % Add text label above bar
    text(1.1,39.81,'20 km','HorizontalAlignment','center','FontSize',16)
    text(lon0+deg_lon/2,lat0,'Mallorca','HorizontalAlignment','center','FontSize',16)

    hold off

    legend([p1, p2], {'Glider', 'CTD'},'Location','southwest')
    RT = load("rt_colormaps.mat");
    colormap(ax1, RT.rt_colormaps.section2);

return

clf
%% Model
% nexttile(2)
% ax = axes('Position',[0.58 0.08 0.40 0.35]);
ax4 = axes('Position',[0.54 0.07 0.40 0.48]);
hour = 6;

    WMOP_path_v6 = '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Simulations/WMOP_REANALYSIS_FaSt-SWOT_2023_v9/';

    dir_WMOP_v6 = dir(fullfile([WMOP_path_v6  '*.nc']));
    [~, data_v6] = read_nc_file_struct([WMOP_path_v6 dir_WMOP_v6(16).name]);
    data_v6.time = datenum(seconds(data_v6.ocean_time) + datenum(1968,5,23));

    % rho  = gsw_rho(data_v6.salt, data_v6.temp, 100);
    rho  = gsw_rho(data_v6.salt, data_v6.temp, 40);

    %% tips for changing density isopycnals,
    %% Use selected density contour levels for the plotted section.

    Var_rho = permute(rho(:,:,:,6), [2 1 3])-1000;       % now size = 170 x 252 x 38

    S = data_v6.salt(:,:,:,6);
    Var_S = permute(S, [2 1 3]);

    lon = double(data_v6.lon);
    lat = double(data_v6.lat);
    depth = -double(data_v6.depth);

    [Lon, Lat, Depth] = meshgrid(lon, lat, depth);   % size = 170 x 252 x 38

    slice(Lon, Lat, Depth, Var_S, [], [39.8], [])
    hold on
    slice(Lon, Lat, Depth, Var_S, [2.2], [], [])
    hold on
    slice(Lon, Lat, Depth, Var_S, [], [], 0)
    shading interp;
    ssha = (data_v6.zeta(:,:,hour) - nanmean(data_v6.zeta(:,:,hour)))';

    hold on
    % contour(data_v6.lon, data_v6.lat, ...
    % (data_v6.zeta(:,:,hour) - nanmean(data_v6.zeta(:,:,hour)))',[-0.3:0.01:0.3],'color',[.5 .5 .5],'LineWidth',1.5)
    hold on
    h1 = contourslice(Lon, Lat, Depth, Var_rho, [], [39.8], [], 28.8:0.05:29.1); %28.8:0.05:30.3);
    set(h1, 'EdgeColor', 'k','LineWidth',1.5);
    hold on
    h2 = contourslice(Lon, Lat, Depth, Var_rho, [2.2], [], [], 28.8:0.05:29.1);
    set(h2, 'EdgeColor', 'k','LineWidth',1.5);

    axis([0.80 2.2 39.8 40.3 -500 0]); grid on
    view([20 35])

    % caxis([28.2 29.4])
    caxis([38.2 38.6])

    % RT = load("rt_colormaps.mat");
    colormap(RT.rt_colormaps.sst2);
    % cmocean('delta')
    box on;
    set(gca,'FontSize',16,'Layer','top')

    hold on

    text(-0.05,1.78,'(c)'...
        ,'Units','normalized','FontWeight','bold','FontSize',16)
    text(-0.05,0.75,'(d)'...
        ,'Units','normalized','FontWeight','bold','FontSize',16)
    title(' Model Salinity (6 May 2023 5:00)'...
        ,'Units','normalized','FontWeight','bold','FontSize',16)

    text(3.2,39.2,'Lat (\circN)','Color','k',...
        'FontSize',18,...
        'Rotation',55)

    text(2.1,38.95,'Lon (\circE)','Color','k',...
        'FontSize',18,...
        'Rotation',-5)

    text(1.9,39.51,'28.95','Color','k',...
        'FontSize',14)

    text(1.9,39.61,'28.9','Color','k',...
        'FontSize',14)

    plot3([0.8 0.8],[39.8 39.8],[-500 0], 'k-','LineWidth',0.5)

cb = colorbar('southoutside');   % horizontal colorbar
cb.Position = [0.4 0.07 0.12 0.012];  % [left bottom width height]
cb.Label.String = 'Salinity (g kg^{-1})';
cb.Label.FontSize = 14;

% close all
% print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig1_left'], '-dpng', '-r300');
% print(['/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig1_right'], '-dpng', '-r300');
