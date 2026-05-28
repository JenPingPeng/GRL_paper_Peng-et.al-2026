# GRL_paper_Peng-et.al-2026
Codes to accompany "Beyond the Surface: Synergistic Assimilation of SWOT and In-Situ Observations to Resolve a Small-scale Eddy"
========================================================================================================================

This repository contains MATLAB scripts used to reproduce the main figures and supporting-information diagnostics for the GRL manuscript. The scripts retain the analysis workflow used for the submitted paper and require the local data paths, MATLAB toolboxes, and helper functions defined within each file.

Main figure scripts
-------------------
- `GRL_Fig1.m`: Generates Figure 1a-d, showing the observational context and the data-assimilation framework.

- `GRL_Fig2_v2.m`: Generates Figure 2a-f, including 2D SSH maps, temporal evolution of the SSH reconstruction, and assimilation-period SSH RMSE diagnostics.

- `GRL_Fig3_v2.m`: Generates Figure 3a-o, including across-eddy SSH, temperature, and salinity sections along the glider transect.

- `GRL_Fig4_SWOT_plot2.m`: Generates Figure 4a-c, showing 2D SSH maps on 9 May 2023 from SWOT DA, Joint DA, and SWOT observations, with ADCP tracks overlaid.

- `GRL_Fig4_ADCP.m`: Generates Figure 4d-f, showing cross-eddy velocity transects from SWOT DA, Joint DA, and independent ADCP observations.

- `GRL_Fig4_SWOT_stats.m`: Generates Figure 4g, showing forecast-period SSH RMSE against SWOT and the associated wind-stress comparison.

- `GRL_Fig4_RMSE_Glider.m`: Computes glider-based RMSE statistics used in the Figure 4 hydrographic validation summary.

- `GRL_Fig4_CTD_RMSE.m`: Computes CTD-based RMSE statistics used in the Figure 4 hydrographic validation summary.

- `GRL_Fig4_histo_glider_CTD.m`: Generates Figure 4h, summarizing glider and CTD RMSE for temperature and salinity.

Supporting-information script
-----------------------------
- `GRL_supporting_cold_waters.m`: Generates the supporting cold-water and salinity-evolution diagnostics at approximately 300 m, used to support the lateral-advection interpretation discussed in the Supporting Information.

Workflow
--------
Run `GRL_Fig1.m`, `GRL_Fig2_v2.m`, and `GRL_Fig3_v2.m` directly to reproduce Figures 1-3.

Figure 4 is assembled from multiple scripts. Run `GRL_Fig4_SWOT_plot2.m` for panels a-c, `GRL_Fig4_ADCP.m` for panels d-f, `GRL_Fig4_SWOT_stats.m` for panel g, and `GRL_Fig4_histo_glider_CTD.m` for panel h. Run `GRL_Fig4_RMSE_Glider.m` and `GRL_Fig4_CTD_RMSE.m` before `GRL_Fig4_histo_glider_CTD.m` to obtain the RMSE values summarized in panel h.

Run `GRL_supporting_cold_waters.m` to reproduce the supporting diagnostics for the temporal evolution of the subsurface cold and fresh anomaly.

Requirements
------------
The scripts require the same local data folders, MATLAB toolboxes, helper functions, and auxiliary files used in the submitted analysis, including TEOS-10, cmocean, the Sea-Bird toolbox, `read_nc_file_struct`, and local coastline and colormap files.

The scripts use local absolute paths, mainly under `/media/jpeng/...`, which are intentionally preserved for reproducibility in the original analysis environment.
