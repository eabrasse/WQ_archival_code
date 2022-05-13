This folder contains the code used by Elizabeth Brasseale to produce the results in the manuscript, "Performance of a one-dimensional model of wave-driven nearshore alongshore tracer transport and decay."

Last updated 4/29/2022

***********************************
Module with miscellaneous tools used throughout
***********************************

"wqfun.py"

***********************************
Code used to make figures
***********************************
"plot_cside_map_paper_B.py" Figure 1
uses
    -grid data '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
    -wave buoy location '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/shoreline_variables_SZ_2017–2018.p'
    -surface dye snapshot '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/surface_dye_01_July.p'

"plot_vel_overview_compare.py" Figure 2
uses
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
    -wave information at wave buoy: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/shoreline_variables_SZ_2017–2018.p'

"plot_alongshore_dye_stats_C.py" Figure 3
uses
    -1D model output from ['CSIDE_2017_2018','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd','AV_recycled_tuned_C0']
	(note: CSIDE_2017_2018 not a 1D model, but COAWST output reformatted for comparison w/ 1D models)

"plot_July2017_E.py" Figure 4
uses
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
    -1D model output from ['U_isobath5m_sawc_autotune_kd']
	
"plot_compare_dye_paper_B.py" Figure 5
uses
    -1D model output from ['CSIDE_2017_2018','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd']
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'

"plot_compare_corrs_paper.py" Figure 6
note: metrics are calculated online in the plotting code and not saved
uses
    -1D model output from ['CSIDE_2017_2018','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd','AV_recycled_tuned_C0']
	
"plot_dyehist_at_IBP_D.py" Figure 7
uses
    -1D model output from ['U_isobath5m_sawc_autotune_kd']
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
	
"plot_dye_at_IBP.py" Figure 8
uses
    -1D model output from ['U_isobath5m_sawc_autotune_kd']
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'

"plot_dye_difference_paper_B.py" Figure 9
uses
    -1D model output from ['CSIDE_2017_2018','U_isobath5m_sawc_autotune_kd']
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
	
"plot_dye_sampling_paper.py" Figure 10
uses
    -1D model output from ['CSIDE_2017_2018','U_isobath5m_sawc_autotune_kd']
    -COAWST model extracted and processed data: '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'


***********************************
Data extraction from COAWST model
***********************************

"x_cside_vars_surfzone_calc.py" Extract wave properties at offshore buoy location for 2017 and 2018 (I only used 2017) Also extracts dye and velocity to estimated surf zone width, but I ultimately didn't use these extractions

"x_cside_vars_calc.py" Extract dye and velocity to 5m isobath for 2017 and 2018 wave properties at offshore buoy location. I only used dye and velocity for 2017 and used "x_cside_vars_surfzone_calc.py" for wave properties.

"x_surface_dye_July.py" Extract snapshots during dye plume for example. Ultimately only used July 11 snapshot for Fig 1b.


***********************************
COAWST data processing
***********************************

"check_rshore.py" Rotate and interpolate velocities

"plot_10m_shorenormal.py" Identifies shoreline indices so it doesn't have to be done online

***********************************
1D model
***********************************

"shoreline_model_tune_PB_in.py" iterative scheme using the same base as "shoreline_model.py". Start with three guesses and it plots the WSS of the three models north of PB. It suggests a next guess based on which results from the initial guesses were closest, but you can override with a value of your own. Once WSS improvement is acceptable, you save output and the figure.

"shoreline_model_cside_tune_PB_in.py" same iterative scheme as "shoreline_tune_PB_in.py" but uses the alongshore-varying base code from "shoreline_model_cside_input.py"

"shoreline_model.py" solves upwind dye advection scheme + dye loss using uniform alongshore velocities derived from wave properties at offshore buoy
Uses
    - wave properties from 'shoreline_variables_SZ_2017–2018.p'
    - grid size from 'shoreline_dye_waves_05m_interp.p' (why different? I can't remember...)
    - kd derived from COAWST results in "plot_alongshore_dye_stats_C.py"
    - PB_in derived from "shoreline_model_tune_PB_in."

"shoreline_model_cside_input.py" solves upwind dye advection scheme + dye loss using alongshore-varying alongshore velocities extracted from COAWST  and processed (rotated and interpolated)
Uses
    - processed velocity and shape of dye array extracted from COAWST from "shoreline_dye_waves_uv_05m_interp_sm100_SA.p"

