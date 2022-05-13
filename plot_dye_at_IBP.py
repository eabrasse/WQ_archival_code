#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot time series of nearshore dye from COAWST and 1D model at Imperial Beach
Demonstrating when dye exceeds cutoff of 5e-4 and identification of four conditions
comparing binary results from the two models: True Positive, False Positive, False Negative,
and True Negative.

The 1D model used here is this one:
1D model
    with velocity calculated at 5m isobath using SAWC approx (see wqfun.wavebuoy_to_5m_isobath(Dbuoy))
    kd = 1.3e-5
    PB_in = 0.008
"""

# setup
import os
import sys
import wqfun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean as cmo
import pickle
import numpy as np

from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import skill_metrics as sm
import matplotlib.colors
import matplotlib.patheffects as PathEffects

# load in color maps
c10 = plt.get_cmap('tab10',10)

TN = 'white'
FP = 'salmon'
TP = 'rebeccapurple'
FN = 'lightskyblue'

four_color_cmap = matplotlib.colors.ListedColormap([TN,FP,FN,TP])
    
plt.close('all')

# load in 1D model data
home = '/Users/elizabethbrasseale/Projects/Water quality/'
model_name_list = ['U_isobath5m_sawc_autotune_kd']
D = wqfun.get_shoreline_models(model_name_list)

# get location of Imperial beach
beach_name_list = ['IB']
beach_list = wqfun.get_beach_location(beach_name_list)

# find dye concentrations from 1D model at Imperial Beach
yi = np.argmin(np.abs(D['U_isobath5m_sawc_autotune_kd']['y']-beach_list['IB']['r']))
dye_1D = D['U_isobath5m_sawc_autotune_kd']['dye'][:,yi]

# load in COAWST model data
cside_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Dcside = pickle.load(open(cside_fn,'rb'))
rshore_old = Dcside['rshore_old'][:]
roi = np.argmin(np.abs(rshore_old-beach_list['IB']['r']))
dye_COAWST = Dcside['dye_01'][:,roi]

# convert time series from 'seconds since Jan 1, 1999' to datetimes
ot=Dcside['ot'][:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

# define hazard threshold
hazard = 5e-4


# generate figure and axis handles
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,fh))
ax = fig.gca()

lw = 1.2

# plot time series of dye from both models
ax.plot(dt_list,dye_1D,color=c10(0),lw=lw,label='1D model')
ax.plot(dt_list,dye_COAWST,color='k',lw=lw,linestyle='solid',label='COAWST')
ax.set_ylabel(r'Dye concentration')
ax.set_yscale('log')
ax.axhline(y=hazard,color='red',linestyle='dashed') # add cut off line
ax.grid(alpha=1,color='lightgray')
ax.text(0.1,0.7,'SD Bight model',color='k',transform=ax.transAxes)
ax.text(0.1,0.6,'1D model',color=c10(0),transform=ax.transAxes)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d"))
plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax.set_xlabel('Date of 2017',fontweight='normal')

# now calculate where dye > hazard in both models
coawst_boolean = dye_COAWST>hazard
coawst_array = np.array(coawst_boolean,dtype=float)*2 #can either be 0 or 2
dataset_boolean = dye_1D>hazard
dataset_array = np.array(dataset_boolean,dtype=float) #can either be 0 or 1
dye_four_option = coawst_array+dataset_array # can be 0 (both no), 1 (model yes, cside no), 2 (model no, cside yes), or 3 (both yes)
# draw bars at the top of the plot indicating duration of each condition
yinds = np.array([7.5e-3,1e-1])
dye_four_option_stag = 0.5*(dye_four_option[1:]+dye_four_option[:-1])
dye_four_option_stag = np.reshape(dye_four_option_stag,(1,len(dye_four_option_stag)))
pv=ax.pcolormesh(dt_list,yinds,dye_four_option_stag,cmap=four_color_cmap,shading='flat')
ax.set_ylim([1e-4,1e-2])
ax.set_xlim([datetime(2017,6,1),datetime(2017,9,1)])

# make a legend to interpret colors
fs = 10
h = 1.1
ax.text(0.98,0.9,'False negative',color=FN,transform=ax.transAxes,ha='right',va='top',fontsize=fs,fontweight='bold')
ax.text(0.98,0.84,'False positive',color=FP,transform=ax.transAxes,ha='right',va='top',fontsize=fs,fontweight='bold')
ax.text(0.98,0.78,'True positive',color=TP,transform=ax.transAxes,ha='right',va='top',fontsize=fs,fontweight='bold')
txt_TN = ax.text(0.98,0.72,'True negative',color=TN,transform=ax.transAxes,ha='right',va='top',fontsize=fs,fontweight='bold')
txt_TN.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])

# show or save plot
plt.tight_layout()
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_08_C.eps'
plt.savefig(outfn,format='eps')
