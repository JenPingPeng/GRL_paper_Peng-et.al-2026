# GRL_paper_Peng-et.al-2016
Codes to accompany &quot;Synergistic Assimilation of SWOT and In-situ Observations for the 3D Representation of a Small-Scale Intrathermocline Eddy&quot;

GRL_Fig1.m → Computations and generation of plots for Figure 1 (panels a-d): observational context and data assimilation framework.

GRL_Fig2.m → Computations and generation of plots for Figure 2 (panels a-p): 2D SSH maps, temporal evolution, and RMSE diagnostics.

GRL_Fig3.m → Computations and generation of plots for Figure 3 (panels a-o): across-eddy SSH, temperature, and salinity sections along the glider transect.

GRL_Fig4_SWOT_plots_v2.m → Computations and generation of plots for Figure 4 (panels a-c): 2D SSH maps on 9 May, including ADCP tracks.

GRL_Fig4_ADCP.m → Computations and generation of plots for Figure 4 (panels d-f): cross-eddy velocity transects from SWOT DA, Joint DA, and ADCP.

GRL_Fig4_SWOT_stats.m → Computations and generation of plots for Figure 4g: forecast-period SSH RMSE time series against SWOT, with wind stress comparison.

GRL_Fig4_RMSE_Glider.m → Computes glider RMSE used for the Figure 4 hydrographic
summary.

GRL_Fig4_RMSE_CTD.m → Computes CTD RMSE used for the Figure 4 hydrographic summary.

GRL_Fig4_histo_glider_CTD.m → Computations and generation of plots for Figure 4h: bar chart summarizing glider and CTD RMSE for temperature and salinity.

Workflow
- Run GRL_Fig1.m, GRL_Fig2.m, and GRL_Fig3.m directly to reproduce Figures 1-3.
  
- Figure 4 is assembled from multiple scripts: GRL_Fig4_SWOT_plots_v2.m for panels a-c, GRL_Fig4_ADCP.m for panels d-f, GRL_Fig4_SWOT_stats.m for panel g, and GRL_Fig4_histo_glider_CTD.m for panel h.

- Run GRL_Fig4_RMSE_Glider.m and GRL_Fig4_RMSE_CTD.m before GRL_Fig4_histo_glider_CTD.m to obtain the RMSE values summarized in panel h.

- The scripts are kept in their submitted analysis form and use the local absolute paths defined inside each MATLAB file.
