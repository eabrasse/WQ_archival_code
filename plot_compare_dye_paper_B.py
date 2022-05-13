#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot year-long nearshore dye concentrations from
COAWST and two 1D nearshore models for comparison

The two 1D nearshore models used are these:
1D model
    with velocity calculated at 5m isobath using SAWC approx (see wqfun.wavebuoy_to_5m_isobath(Dbuoy))
    kd = 1.3e-5
    PB_in = 0.008
1DR model
    with velocity calculated at 5m isobath using Rayleigh friction approx (divide Sxy by linear fit)
                             note: for Rayleigh friction, velocity at 5m isobath same as at wavebuoy
    kd = 1.3e-5
    PB_in = 0.010
"""

# setup
import os
import sys
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
import string

# define a few things for plotting
atoz = string.ascii_lowercase
c10 = plt.get_cmap('tab10',10)
# get beach locations, only really need PB
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)
props = dict(boxstyle='round', fc='white',ec='None',alpha=.5)
    
plt.close('all')



# load in models
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd']

D = wqfun.get_shoreline_models(model_name_list)
CSIDE = D['CSIDE_SAsm100']
dataset_list = [CSIDE]
for model_name in model_name_list[1:]:
    dataset_list.append(D[model_name])

dataset_list[0]['label'] = 'COAWST model'
dataset_list[-1]['label'] = '1D model'
y = CSIDE['y']-beach_list['PB']['r']

# open up data extracted from CSIDE just to use time series
home = '/Users/elizabethbrasseale/Projects/Water quality/'
cside_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Dcside = pickle.load(open(cside_fn,'rb'))

# trim the first few days of the dataset
# data initially 0 everywhere, and logarithmic color scale plots
# it in white as "bad data"
# it's distracting and only happens for the first few days
# another way to solve this problem would be to artificially
# nudge the data by 1e-12 or something else super small, so it's not 0
t0 = 168
t1 = None
ot=Dcside['ot'][:]
dt_list = []
for ott in ot[t0:t1]:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))


# colorscale plotting parameters
dyemax = 5e-2
dyemin = 1e-4
vmax = dyemax
vmin = dyemin
colmap = cmo.cm.matter
normal = matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax)

# generate figure, subplot axes using gridspec, and handles
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,1.5*fh))
nds = len(dataset_list)
gs = GridSpec(nds,1)
axs_var = [] # make handles indexable for flexibility
for i in range(nds):
    ax = fig.add_subplot(gs[i])
    axs_var.append(ax)

#PLOT MODEL DYE
plot_labels = ['SD Bight model', '1D model', '1DR model']
ylim = [0,y.max()*0.001]
count=0
colorcount = 0
for dataset in dataset_list:
    
    ax_vel = axs_var[count] #find axis
    
    # plot dye
    pv = ax_vel.pcolormesh(dt_list,0.001*y,np.transpose(dataset['dye'][t0:t1,:]),norm=normal,cmap=colmap,shading='nearest')
    ax_vel.contour(dt_list,0.001*y,np.transpose(dataset['dye'][t0:t1,:]),levels=[5e-4],linestyles=['dashed'],linewidths=[0.5],colors=['k'])

    # add label
    ax_vel.text(0.01,0.95,atoz[count]+') '+plot_labels[count],transform=ax_vel.transAxes,ha='left',va='top',bbox=props)
    ax_vel.set_ylabel('y (km)')
    ax_vel.set_ylim(ylim)
    wqfun.add_beaches_to_ax_RHS(ax_vel,TJRE_line=False,xaxscale=0.02)
    
    # plot colorbar once
    if count==0:
        
        cbaxes = inset_axes(ax_vel, width="30%", height="6%", loc=6,bbox_transform=ax_vel.transAxes,bbox_to_anchor=(0.01,0.18,1,1))
        cb = fig.colorbar(pv, cax=cbaxes, orientation='horizontal')
        
    
    # plot x-axis on bottom axis only; otherwise, hide x label
    if count==nds-1:
        ax_vel.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d"))
        plt.setp( ax_vel.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
        ax_vel.set_xlabel('Date of 2017',fontweight='normal')
    else:
        ax_vel.set_xticklabels([''])
        ax_vel.get_xaxis().set_visible('False')
        
    count+=1

# show or save plot
plt.subplots_adjust(bottom=0.12,hspace=0.1,top=0.98)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_06_E.jpg'
plt.savefig(outfn,format='jpg',dpi=600)
