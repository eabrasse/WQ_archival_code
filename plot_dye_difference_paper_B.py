#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot where 1D model and 3D model overlap in presence of dye > 5e-4
  True Positive: dye > 5e-4 in both models
  False Positive: dye > 5e-4 only in 1D model
  False Negative: dye > 5e-4 only in 3D model
  True Negative: dye > 5e-4 in neither model
This will be demonstrated as a function of alongshore distance and time
then as a cumulative percentage of each condition (adding to 100%) as a function of distance
then as a percent agreement (both or neither > 5e-4) as a function of distance

The 1D model used here is
1D model
    with velocity calculated at 5m isobath using SAWC approx (see wqfun.wavebuoy_to_5m_isobath(Dbuoy))
    kd = 1.3e-5
    PB_in = 0.008
"""

# setup
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
import matplotlib.patheffects as PathEffects

# useful plotting tools, colormap and test box props
TN = 'white'
FP = 'salmon'
TP = 'rebeccapurple'
FN = 'lightskyblue'

four_color_cmap = matplotlib.colors.ListedColormap([TN,FP,FN,TP])

props = dict(boxstyle='round', facecolor='white')

plt.close('all')

# load in 1D model and 3D model
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd']
# model_name_list = ['CSIDE_SAsm100','AV_recycled_tuned_C0']
D = wqfun.get_shoreline_models(model_name_list)
dataset = D[model_name_list[-1]]
CSIDE = D['CSIDE_SAsm100']
nt,nj = dataset['dye'].shape

# load in COAWST data to build time series
home = '/Users/elizabethbrasseale/Projects/Water quality/'
cside_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Dcside = pickle.load(open(cside_fn,'rb'))
ot=Dcside['ot'][:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

t0 = 0
t1 = None

#load in beach locations
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)
PB = beach_list['PB']
y = CSIDE['y'][:]-PB['r']


# generate figure and subplots
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,1.5*fh))
gs = GridSpec(2,2)

# calculate 4 conditions
hazard = 5e-4
cside_boolean = CSIDE['dye'][t0:t1,:]>hazard
cside_array = np.array(cside_boolean,dtype=float)*2 #can either be 0 or 2
dataset_boolean = dataset['dye'][t0:t1,:]>hazard
dataset_array = np.array(dataset_boolean,dtype=float) #can either be 0 or 1
dye_four_option = cside_array+dataset_array # can be 0 (both no), 1 (model yes, cside no), 2 (model no, cside yes), or 3 (both yes)

# plot 4 conditions as function of distance and time
ax_bin = fig.add_subplot(gs[0,:])
pv=ax_bin.pcolormesh(dt_list[t0:t1],0.001*y,np.transpose(dye_four_option),cmap=four_color_cmap,shading='nearest')
ax_bin.set_ylabel('y (km)')
ax_bin.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d"))
plt.setp( ax_bin.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_bin.set_xlabel('Date of 2017',fontweight='normal')
ax_bin.text(0.05,0.95,'a)',transform=ax_bin.transAxes,ha='left',va='top')

ylim = ax_bin.get_ylim()
ax_bin.set_ylim([0,ylim[1]])

# now estimate percentage of both in agreement
number_correct = np.sum(cside_boolean==dataset_boolean,axis=0)
percent_correct = 100*number_correct/nt

# plot percent of time steps in agreement
ax_pc = fig.add_subplot(gs[1:,1])
ax_pc.plot(percent_correct,0.001*y,color='k',linestyle='solid',lw=2,label='true pos&neg/all')
ax_pc.grid(color='lightgray',which='both')
axtickrange = range(80,101,10)
ax_pc.set_xticks(axtickrange)
ax_pc.set_xticklabels([f'{r}'+r'%' for r in axtickrange])
plt.setp(ax_pc.get_yticklabels(), visible=False)
ax_pc.text(0.05,0.95,'c)',transform=ax_pc.transAxes,ha='left',va='top')
ax_pc.set_xlim([70,100])
ax_pc.set_xlabel('Percent of time steps\nin agreement',ha='center')

#make stacked bar plot of percent of time steps in agreement as a function of distance
ax_sb = fig.add_subplot(gs[1:,0],sharey=ax_pc)

number_true_positive = np.sum((cside_boolean==1)&(dataset_boolean==1),axis=0)
number_true_negative = np.sum((cside_boolean==0)&(dataset_boolean==0),axis=0)
number_false_positive = np.sum((cside_boolean==0)&(dataset_boolean==1),axis=0)
number_false_negative = np.sum((cside_boolean==1)&(dataset_boolean==0),axis=0)

norm = 100/nt
ax_sb.barh(0.001*y,norm*number_true_positive,color=TP)
ax_sb.barh(0.001*y,norm*number_false_positive,left=norm*number_true_positive,color=FP)
ax_sb.barh(0.001*y,norm*number_false_negative,left=norm*(number_true_positive+number_false_positive),color=FN)
ax_sb.barh(0.001*y,norm*number_true_negative,left=norm*(number_true_positive+number_false_positive+number_false_negative),color=TN)
ax_sb.set_ylabel('y (km)')
axtickrange = range(0,101,20)
ax_sb.set_xticks(axtickrange)
ax_sb.set_xticklabels([f'{r}'+r'%' for r in axtickrange])
ax_sb.set_ylim([0,0.001*y.max()])
ax_sb.axhline(y=0.001*(beach_list['TJRE']['r']-beach_list['PB']['r']),color='yellowgreen',linestyle='dashed',linewidth=1.0)
ax_sb.text(0.05,0.95,'b)',transform=ax_sb.transAxes,ha='left',va='top')
ax_sb.set_xlabel('Percent of time steps\nof all conditions',ha='center')

# add legend
fs = 10
ax_sb.text(0.5,0.95,'False Negative',color=FN,transform=ax_sb.transAxes,ha='left',va='top',fontsize=fs,fontweight='bold')
ax_sb.text(0.5,0.87,'False Positive',color=FP,transform=ax_sb.transAxes,ha='left',va='top',fontsize=fs,fontweight='bold')
ax_sb.text(0.5,0.79,'True Positive',color=TP,transform=ax_sb.transAxes,ha='left',va='top',fontsize=fs,fontweight='bold')
txt_TN = ax_sb.text(0.5,0.71,'True Negative',color=TN,transform=ax_sb.transAxes,ha='left',va='top',fontsize=fs,fontweight='bold')
txt_TN.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])

# add beach locations
wqfun.add_beaches_to_ax_RHS(ax_pc,TJRE_line=True,xaxscale=0.05)
wqfun.add_beaches_to_ax_RHS(ax_bin,TJRE_line=True,xaxscale=0.02)
wqfun.add_beaches_to_ax_RHS(ax_sb,TJRE_line=True,xaxscale=0.02*10/3)

# show or save plot
plt.tight_layout()
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_09_D.jpg'
plt.savefig(outfn,format='jpg',dpi=600)
