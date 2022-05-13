#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a single dye plume event in the COAWST model alongside
the nearshore alongshore velocities and dye in the 1D model

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

#define some formatting parameters right away
props = dict(boxstyle='round', facecolor='white',edgecolor='white')
c10 = plt.get_cmap('tab20c',10)

fs_small = 10
fs_big = 12
    
plt.close('all')

# load in files
home = '/Users/elizabethbrasseale/Projects/Water quality/'
wave_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p' # has dye & velocity from COAWST
Dw = pickle.load(open(wave_fn,'rb'))

#convert "seconds since Jan 1, 1999" to list of datetimes for plotting
ot = Dw['ot'][:]
dt_list0 = []
for ott in ot:
    dt_list0.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

# find which indexes in list to use for limiting time period for plotting
# it's more legible to do this with datetimes, but finding the nearest ind
# in a list of non-numeric data is a little annoying
dt0 = datetime(2017,7,8,12)
dt1 = datetime(2017,7,12)
td_list0 = [np.abs(dt-dt0) for dt in dt_list0]
t0 = np.argmin(td_list0)
td_list1 = [np.abs(dt-dt1) for dt in dt_list0]
t1 = np.argmin(td_list1)
dt_list=dt_list0[t0:t1]

# Note: I'm using the NON-interpolated dye and v,
# since interpolation introduces errors
# I use the interpolated data when doing direct
# comparisons with the 1D model output, but that's not
# necessary here.
rshore = Dw['rshore_old'][:] # rshore_old has irregular grid spacing, rshore is interpolated onto regular grid
Vco = np.mean(Dw['v_rot'][t0:t1,:],axis=1)

# load in 1D model
model_name_list = ['U_isobath5m_sawc_autotune_kd']
D = wqfun.get_shoreline_models(model_name_list)
model = D[model_name_list[0]]
Vmod = model['v'][t0:t1,0]

# count indexes in both space and time
nt,nj = Dw['dye_01'][t0:t1,:].shape

# initialize figure and axes
fw,fh = wqfun.gen_plot_props(fs_small=fs_small,fs_big=fs_big)
fig = plt.figure(figsize=(fw*2,1.5*fh))
gs = GridSpec(8,1)
# note: I name figure handles as I go in the subsequent code chunks

# load in beach information to convert y into 'km from PB'
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)
PB = beach_list['PB']

# plot time series of velocities
ax_sxy = fig.add_subplot(gs[:2,0])
ax_sxy.plot(dt_list,Vco,color='k',label=r'$\bar{v}$')
ax_sxy.plot(dt_list,Vmod,color=c10(0),label=r'$v^{*}$')
ax_sxy.text(0.1,0.9,r'$\bar{v}_{C}$',color='k',transform=ax_sxy.transAxes,va='top',ha='left',bbox=props,fontsize=fs_small)
ax_sxy.text(0.1,0.7,r'$v_{1D}$',color=c10(0),transform=ax_sxy.transAxes,va='top',ha='left',bbox=props,fontsize=fs_small)
ax_sxy.axhline(y=0,color='k',linewidth=0.5)
ax_sxy.get_xaxis().set_visible(False)
ax_sxy.set_xlabel('')
ax_sxy.text(0.025,0.9,'b)',transform=ax_sxy.transAxes,fontsize=fs_big,ha='left',va='top')
ax_sxy.set_ylabel('Velocity\n(m/s)')
ax_sxy.set_xlim([dt_list[0],dt_list[-1]])

# plot COAWST dye
ax_dye = fig.add_subplot(gs[2:5,0])#,sharex=ax_sxy)
p3 = ax_dye.pcolormesh(dt_list,0.001*(rshore-PB['r']),np.transpose(Dw['dye_01'][t0:t1,:]),norm=LogNorm(vmin=1e-4,vmax=5e-2),cmap=cmo.cm.matter,shading='nearest')
ax_dye.contour(dt_list,0.001*(rshore-PB['r']),np.transpose(Dw['dye_01'][t0:t1,:]),levels=[5e-4],linestyles=['dashed'],linewidths=[0.5],colors=['k'])
ax_dye.get_xaxis().set_visible(False)
ax_dye.set_xlabel('Time',fontweight='normal')
ax_dye.set_ylabel('y (km)',fontweight='normal')
ax_dye.text(0.025,0.9,'a)',transform=ax_dye.transAxes,fontsize=fs_big,ha='left',va='top')
ax_dye.set_xlim([dt_list[0],dt_list[-1]])
ylim = ax_dye.get_ylim()
ax_dye.set_ylim([0,ylim[1]])
wqfun.add_beaches_to_ax_RHS(ax_dye,TJRE_line=False,xaxscale=0.02)

# plot 1D model dye
ax_model = fig.add_subplot(gs[5:,0])#,sharex=ax_sxy)
p3 = ax_model.pcolormesh(dt_list,0.001*(model['y']-PB['r']),np.transpose(model['dye'][t0:t1,:]),norm=LogNorm(vmin=1e-4,vmax=5e-2),cmap=cmo.cm.matter,shading='nearest')
ax_model.contour(dt_list,0.001*(model['y']-PB['r']),np.transpose(model['dye'][t0:t1,:]),levels=[5e-4],linestyles=['dashed'],linewidths=[0.5],colors=['k'])
ax_model.set_xticks([datetime(2017,7,7,12),datetime(2017,7,8,12),datetime(2017,7,9,12),datetime(2017,7,10,12),datetime(2017,7,11,12)])
ax_model.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d"))
plt.setp( ax_model.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_model.set_xlabel('Date of 2017',fontweight='normal')
ax_model.set_ylabel('y (km)',fontweight='normal')
ax_model.text(0.025,0.9,'c)',transform=ax_model.transAxes,fontsize=fs_big,ha='left',va='top')
ax_model.set_xlim([dt_list[0],dt_list[-1]])
ylim = ax_model.get_ylim()
ax_model.set_ylim([0,ylim[1]])
cbaxes = inset_axes(ax_dye, width="40%", height="8%", loc=6,bbox_transform=ax_dye.transAxes,bbox_to_anchor=(0.1,0.32,1,1))
cb = fig.colorbar(p3, cax=cbaxes, orientation='horizontal')
cbaxes.tick_params(axis='both', which='major', labelsize=9)
wqfun.add_beaches_to_ax_RHS(ax_model,TJRE_line=False,xaxscale=0.02)

# add indicator on all axes of when ('where' along x-axis) snapshot in Fig 1b occurs
for ax in ax_sxy,ax_dye,ax_model:
    ax.axvline(x=datetime(2017,7,11,12),color='k',linestyle='dashed')

# show or save figure
plt.subplots_adjust(bottom=0.13,left=0.13,hspace=0.6,top=0.98)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_05_F.jpg'
plt.savefig(outfn,format='jpg',dpi=600)
