#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot comparison of predicted daily beach advisory conditions in the 3D model
with daily beach advisory conditions in the 1D model and simulated weekly sampling
The 1D model used is this one:
1D model
    with velocity calculated at 5m isobath using SAWC approx (see wqfun.wavebuoy_to_5m_isobath(Dbuoy))
    kd = 1.3e-5
    PB_in = 0.008
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
import netCDF4 as nc
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import wqfun
import matplotlib.colors

# useful plotting tools
c10 = plt.get_cmap('tab10',10)
props = dict(boxstyle='round', facecolor='white')
plt.close('all')

# load in 1D models
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd']
D = wqfun.get_shoreline_models(model_name_list)
dataset = D[model_name_list[-1]]
CSIDE = D['CSIDE_SAsm100']

# load in time series
home = '/Users/elizabethbrasseale/Projects/Water quality/'
cside_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Dcside = pickle.load(open(cside_fn,'rb'))
nt,nj = dataset['dye'].shape
t0 = 0
t1 = None
ot=Dcside['ot'][:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))
    
# load in PB location and convert y to 'km from PB'
beach_name_list = ['PB']
beach_list = wqfun.get_beach_location(beach_name_list)
PB = beach_list['PB']
y = CSIDE['y'][:]-PB['r']

# generate figure
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(fw,fh))

# calculate sampling accuracy
percent_correct_sampling = np.zeros((nj))
percent_correct_model = np.zeros((nj))
hazard = 5e-4
for j in range(nj):
    beach_open_sampling = np.ones((nt))
    # sample once a week
    for t in range(0,nt,24*7):
        #default is YES beach is open
        beach_open = 1
    
        #if dye conc is high at the time of sampling
        if CSIDE['dye'][t,j]>hazard:
            # NO beach should not be open
            beach_open = 0

        # if bad water quality, close beach for 48 hours
        # beach_open_sampling[t:(t+48)] = beach_open
        # then reopen until next week
        # beach_open_sampling[(t+48):(t+24*7)] = 1
    
        #or try just closing for the week
        #include delay of processing results
        beach_open_sampling[t+24:(t+8*24)] = beach_open

    beach_open_model = np.ones((nt))
    for t in range(0,nt,24):
        #default is YES beach is open
        beach_open = 1
        # if model says dye concentrations are high at any point that day...
        if np.any(dataset['dye'][t:(t+24),j]>hazard):
            # NO beach is not open
            beach_open = 0
        # beach status will be kept for the whole day
        beach_open_model[t:(t+24)] = beach_open
    
    beach_open_perfect = np.ones((nt))
    for t in range(0,nt,24):
        #default is YES beach is open
        beach_open = 1
        # if model says dye concentrations are high at any point that day...
        if np.any(CSIDE['dye'][t:(t+24),j]>hazard):
            # NO beach is not open
            beach_open = 0
        # beach status will be kept for the whole day
        beach_open_perfect[t:(t+24)] = beach_open
        
    percent_correct_sampling[j] = 100*np.sum(beach_open_sampling==beach_open_perfect)/nt
    percent_correct_model[j] = 100*np.sum(beach_open_model==beach_open_perfect)/nt


# plot percent correct for 1D model and simulated weekly sampling
ax_pc = fig.gca()
ax_pc.plot(percent_correct_model,0.001*y,color=c10(0),label='1D model')
ax_pc.plot(percent_correct_sampling,0.001*y,color=c10(6),label='weekly sampling')
ax_pc.set_xlabel('Percent agreement')
ax_pc.set_xticks([60,70,80,90,100])
ax_pc.set_xticklabels(['60%','70%','80%','90%','100%'])
ax_pc.text(90,10,'1D model forecast',color=c10(0),rotation=60,ha='center',va='center')
ax_pc.text(75,17,'Weekly sampling',color=c10(6),rotation=45,ha='center',va='center')
ax_pc.set_ylabel('y (km)')
plt.setp( ax_pc.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_pc.grid(alpha=1.0,color='lightgray')
# ylim = ax_pc.get_ylim()
# ax_pc.set_ylim([0.1,29])

wqfun.add_beaches_to_ax_RHS(ax_pc,TJRE_line=True,xaxscale=0.05)

# show or save figure
plt.subplots_adjust(bottom=0.2,left=0.18,right=0.85)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_10_D.eps'
plt.savefig(outfn,format='eps')