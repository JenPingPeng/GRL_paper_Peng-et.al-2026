% ============================================================
% Figure 4 release script
% Panels: Figure 4h
% Purpose: Bar chart summarizing glider and CTD RMSE for temperature and salinity.
%
% Notes for release:
% - Run this script for panel (h).
% - The values are entered manually at the top of the script.
% - This script is the final summary panel of Figure 4.
% ============================================================

clear; clc; close all;

%% =========================================================
% USER DATA
%% =========================================================

% Glider mean RMSE
glider_T = [0.47, 0.39];   % [SWOT DA, Joint DA]
glider_S = [0.08, 0.07];

% CTD mean RMSE
ctd_T = [0.30, 0.25];      % [SWOT DA, Joint DA]
ctd_S = [0.11, 0.06];

% Colors
col_swot  = [0.36 0.78 0.33];   % green
col_joint = [0.96 0.56 0.20];   % orange

%% =========================================================
% FIGURE
%% =========================================================
fig = figure('Color','w','Position',[150 120 360 360]);
ax = axes('Position',[0.18 0.12 0.80 0.78]);
hold(ax,'on')

%% =========================================================
% X positions
%% =========================================================
x_gT = 1.0;
x_gS = 2.3;
x_cT = 3.7;
x_cS = 5;

dx    = 0.22;
bar_w = 0.36;

%% =========================================================
% Bars
%% =========================================================
h1 = bar(ax, x_gT-dx, glider_T(1), bar_w, ...
    'FaceColor', col_swot, 'EdgeColor', 'k', 'LineWidth', 0.7);
h1.FaceAlpha = 0.85;

h2 = bar(ax, x_gT+dx, glider_T(2), bar_w, ...
    'FaceColor', col_joint, 'EdgeColor', 'k', 'LineWidth', 0.7);
h2.FaceAlpha = 0.85;

h3 = bar(ax, x_gS-dx, glider_S(1), bar_w, ...
    'FaceColor', col_swot, 'EdgeColor', 'k', 'LineWidth', 0.7);
h3.FaceAlpha = 0.85;

h4 = bar(ax, x_gS+dx, glider_S(2), bar_w, ...
    'FaceColor', col_joint, 'EdgeColor', 'k', 'LineWidth', 0.7);
h4.FaceAlpha = 0.85;

h5 = bar(ax, x_cT-dx, ctd_T(1), bar_w, ...
    'FaceColor', col_swot, 'EdgeColor', 'k', 'LineWidth', 0.7);
h5.FaceAlpha = 0.85;

h6 = bar(ax, x_cT+dx, ctd_T(2), bar_w, ...
    'FaceColor', col_joint, 'EdgeColor', 'k', 'LineWidth', 0.7);
h6.FaceAlpha = 0.85;

h7 = bar(ax, x_cS-dx, ctd_S(1), bar_w, ...
    'FaceColor', col_swot, 'EdgeColor', 'k', 'LineWidth', 0.7);
h7.FaceAlpha = 0.85;

h8 = bar(ax, x_cS+dx, ctd_S(2), bar_w, ...
    'FaceColor', col_joint, 'EdgeColor', 'k', 'LineWidth', 0.7);
h8.FaceAlpha = 0.85;

%% =========================================================
% Axis style
%% =========================================================
set(ax, 'Box','off', ...
        'TickDir','out', ...
        'LineWidth',1.0, ...
        'FontSize',12, ...
        'FontName','Helvetica', ...
        'Layer','top')

xlim([0.4 5.55])
ylim([0 0.60])

ylabel('RMSE (°C, g kg^{-1})', 'FontSize', 13)

set(ax, 'XTick', [x_gT x_gS x_cT x_cS], ...
        'XTickLabel', {'Temp.','Sal.','Temp.','Sal.'})

ax.YGrid = 'on';
ax.GridAlpha = 0.15;

%% =========================================================
% Group labels
%% =========================================================
text(mean([x_gT x_gS]), 0.535, 'Glider', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',14, 'FontWeight','bold');

text(mean([x_cT x_cS]), 0.535, 'CTD', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',14, 'FontWeight','bold');

plot(ax, [x_gT-0.38, x_gS+0.38], [0.515 0.515], 'k-', 'LineWidth', 1.0);
plot(ax, [x_cT-0.38, x_cS+0.38], [0.515 0.515], 'k-', 'LineWidth', 1.0);

%% =========================================================
% Divider
%% =========================================================
plot(ax, [2.95 2.95], [0 0.50], ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);

%% =========================================================
% Legend
%% =========================================================
lgd = legend(ax, [h1 h2], {'SWOT DA','Joint DA'}, ...
    'Location','none', 'Box','off', 'FontSize',11);
lgd.Position = [0.73 0.66 0.18 0.10];

%% =========================================================
% Panel label
%% =========================================================
text(ax, 0.02, 0.98, '(h)', 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================
% Export
%% =========================================================
set(gcf,'Renderer','painters');
% print(fig, '/media/jpeng/Extreme SSD/01_IMEDEA_2025/Figures/GRL/Fig4_hist', '-dpng', '-r300');

