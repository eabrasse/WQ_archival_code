#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code to plot time-averaged dye concentrations from
COAWST model with fit
plus a number of additional models
as of 4/26/2022, it plots
1D model
    with velocity calculated at 5m isobath using SAWC approx (see wqfun.wavebuoy_to_5m_isobath(Dbuoy))
    kd = 1.3e-5
    PB_in = 0.008
1DR model
    with velocity calculated at 5m isobath using Rayleigh friction approx (divide Sxy by linear fit)
                             note: for Rayleigh friction, velocity at 5m isobath same as at wavebuoy
    kd = 1.3e-5
    PB_in = 0.010
1DC model
    with alongshore-varying velocity extracted, rotated, and interpolated from COAWST model
    kd = 1.3e-5
    PB_in = 0.011
"""

# setup
import os
import sys
alp = os.path.abspath('/Users/elizabethbrasseale/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean as cmo
import pickle
import numpy as np
import wqfun
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import skill_metrics as sm

# generate colors from colormap
c20c = plt.get_cmap('tab20c',20)
cC = 'k'
c1D = c20c(0)
c1DC = c20c(6)
c1DR = c20c(9)
c_list = [c1D,c1DR,c1DC]
    
plt.close('all')

props = dict(boxstyle='round', facecolor='white',edgecolor='None')
# load in 1D models and 1Dized COAWST data
#
# Note: CSIDE_SAsm100 is the same data as in file
# home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p',
# but reformatted so that it can be handled like the 1D models.
# see entry for CSIDE_SAsm100 in wqfun.get_shoreline_models()
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd','AV_recycled_tuned_C0']
D = wqfun.get_shoreline_models(model_name_list)
CSIDE = D['CSIDE_SAsm100']

home = '/Users/elizabethbrasseale/Projects/Water quality/'

# load in beach location information to adjust y-axis to 'km from PB'
beach_name_list = ['PB','TJRE']
beach_list = wqfun.get_beach_location(beach_name_list)

y = CSIDE['y']-beach_list['PB']['r']
# data only decays linearly > 5km from PB
# need to identify this index when running the fit
y_ind = np.argmin(np.abs(y-5000))

#generate figure and axis handles
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(1.2*fw,fh))
ax_avg = fig.gca()

# plot time-averaged dye concentrations for all models
color_count = 0
for model_name in model_name_list:
    data = D[model_name]
    if model_name=='CSIDE_SAsm100':
        color='k'
        lw = 3
    else:
        color = c_list[color_count]
        lw = 2
        color_count+=1
    data_mean = np.mean(data['dye'],axis=0)
    
    ax_avg.semilogx(data_mean,0.001*y,color=color,linestyle='solid',linewidth=lw)


# add axis labels and legends
ax_avg.set_xlabel('Time-averaged nearshore dye')
ax_avg.set_ylabel('y (km)')
ax_avg.text(0.4,0.83,r'$\langle C_{\mathrm{C}}\rangle$',color=cC,transform=ax_avg.transAxes,ha='left',fontweight='bold',bbox=props)
ax_avg.text(0.4,0.76,r'$\langle C_{\mathrm{1D}}\rangle$',color=c1D,transform=ax_avg.transAxes,ha='left',bbox=props)
ax_avg.text(0.4,0.69,r'$\langle C_{\mathrm{1DR}}\rangle$',color=c1DR,transform=ax_avg.transAxes,ha='left',bbox=props)
ax_avg.text(0.4,0.62,r'$\langle C_{\mathrm{1DC}}\rangle$',color=c1DC,transform=ax_avg.transAxes,ha='left',bbox=props)

# format axis
ax_avg.grid(color='lightgray',alpha=1.0)
# ax_avg.set_ylim([1,0.001*y[-1]])
ax_avg.set_xlim([5e-5,5e-3])

# calculate & plot linear regression of log of time-averaged COAWST dye
dye = np.mean(CSIDE['dye'],axis=0)
dyelog = np.log(dye[y_ind:]) # if D(y) = D0 * exp(-k_d * y), log(D) = log(D0) - k_d * y
p = np.polyfit(y[y_ind:],dyelog,1) # should be log(D) = p[0] * y + p[1]; p[1] = log(D0), p[0] = -k_d
D0 = np.exp(p[1])
D = np.exp(p[1] + p[0]*y[y_ind:])
ax_avg.semilogx(D,0.001*y[y_ind:],color=cC,linestyle='dashed',lw=1)
L = -1.0 / p[0]

# add beach locations to RHS for reference    
wqfun.add_beaches_to_ax_RHS(ax_avg,TJRE_line=True,xaxscale=0.2)


# show or save plot
plt.tight_layout()
# plt.show(block=False)
# plt.pause(0.1)
# outfn = home+'WQ_plots/paper figures/Figure_04_B.jpg'
# plt.savefig(outfn,format='jpg',dpi=600)
outfn = home+'WQ_plots/paper figures/Figure_04_B.eps'
plt.savefig(outfn,format='eps')