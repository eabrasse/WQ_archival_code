#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot comparing three performance metrics (R, NRMSE, and WSS)
for three models:
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

# load in some useful plotting things (colormap and a-to-z string)
atoz = string.ascii_lowercase

c20c = plt.get_cmap('tab20c',20)
c1D = c20c(0)
c1DC = c20c(6)
c1DR = c20c(9)
c_list = [c1D,c1DR,c1DC]

props = dict(boxstyle='round', facecolor='white',edgecolor='none')

# define performance metrics as functions so they can be called the same way
# and output the same way
def xcorr_solve(mod,obs,nj):
    xcorr = np.zeros((nj))
    for j in range(nj):
        xcorr[j] = np.corrcoef(obs[:,j],mod[:,j])[0,1]
    return xcorr

def rmsd_solve(mod,obs,nj):
    rmsd = np.zeros((nj))
    for j in range(nj):
        rmsd[j] = sm.rmsd(mod[:,j],obs[:,j])
    return rmsd

def nrmsd_solve(mod,obs,nj):
    nrmsd = np.zeros((nj))
    for j in range(nj):
        # obs_mean = np.mean(obs[:,j])
        obs_rms = np.sqrt(np.mean(obs[:,j]**2))
        rmsd = sm.rmsd(mod[:,j],obs[:,j])
        nrmsd[j] = rmsd/obs_rms
    return nrmsd

def wss_solve(mod,obs,nj):
    Wskillscore = np.zeros((nj)) 
    for j in range(nj):
        Wskillscore[j] = wqfun.willmott(mod[:,j],obs[:,j])
    return Wskillscore
    
plt.close('all')


# load in models
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd','U_isobath5m_R_autotune_kd','AV_recycled_tuned_C0']
D = wqfun.get_shoreline_models(model_name_list)
CSIDE = D['CSIDE_SAsm100']
model_list = []
for model_name in model_name_list[1:]:
    locals()[model_name]=D[model_name]
    model_list.append(locals()[model_name])
    
# adjust y axis to 'km from PB'
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)
y = CSIDE['y']-beach_list['PB']['r']
PB_ind = np.argmin(np.abs(y))
ylim = [2,y.max()*0.001]

# define metrics as dicts so that they can be plotted iteratively
# or if you want to plot some smaller selection of them, it's easy
xcorr = {}
xcorr['func'] = xcorr_solve
xcorr['title'] = 'R'
xcorr['xlim'] = [0,1]
xcorr['xticks'] = [0,0.5,1.0]
xcorr['xminorticks'] = [0,0.25,0.5,0.75,1.0]

RMSE = {}
RMSE['func'] = rmsd_solve
RMSE['title'] = 'RMSE'
RMSE['xlim'] = [1e-5,1e-1]
RMSE['xticks'] = [1e-5,1e-3,1e-1]
RMSE['xminorticks'] = [1e-5,1e-4,1e-3,1e-2,1e-1]

NRMSE = {}
NRMSE['func'] = nrmsd_solve
NRMSE['title'] = 'NRMSE'
NRMSE['xlim'] = [0,2]
NRMSE['xticks'] = [0,1,2]
NRMSE['xminorticks'] = NRMSE['xticks']

WSS = {}
WSS['func'] = wss_solve
WSS['title'] = 'WSS'
WSS['xlim'] = [0,1]
WSS['xticks'] = [0,0.5,1.0]
WSS['xminorticks'] = [0,0.25,0.5,0.75,1.0]

# itemize which metrics you want to plot
model_eval_list = [xcorr,NRMSE,WSS]

# get axis shape
nt, nj = CSIDE['dye'].shape
t0 = 0
t1 = nt

# generate figure and axis 
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,fh))
nme = len(model_eval_list)
gs = GridSpec(1,nme)

# include some extra formatting details for 
# different models, if you want to format them differently
#
# warning: lazy coding ahead! if you change the number of models
# you plot, errors will appear here

# 1D model
model_list[0]['col'] = c1D
model_list[0]['lw'] = 1.5
model_list[0]['ls'] = 'solid'
# 1DR model
model_list[1]['col'] = c1DR
model_list[1]['lw'] = 1.5
model_list[1]['ls'] = 'solid'
# 1DC model
model_list[2]['col'] = c1DC
model_list[2]['lw'] = 1.5
model_list[2]['ls'] = 'solid'

# I coded this awkwardly, so that the actual performance metric data I generate is plotted but not easily accessible
# When I wanted to find means or maxes of the metrics, I set 'print_stats=True' and it prints it to the terminal
print_stats =True
psc = 0 #count of number of print stats so that you can stagger them
count=0
# loop through performance metrics, calculating the metric for each model and plotting on each axis before moving
for model_eval in model_eval_list:
    ax = fig.add_subplot(gs[0,count])
    model_eval['ax'] = ax
    for model in model_list:
        model_eval['var'] = model_eval['func'](CSIDE['dye'][t0:t1,:],model['dye'][t0:t1,:],nj)
        model_eval['ax'].plot(model_eval['var'],0.001*y,color=model['col'],lw=model['lw'],linestyle=model['ls'])
        if print_stats:
            
            if model_eval==WSS:
                mean_var = np.mean(model_eval['var'][PB_ind+1:])
                max_var = np.max(model_eval['var'][PB_ind+1:])
                print(model['label'])
                stat_print_string = 'mean '+model_eval['title']+f' = {mean_var:.2f}'+'\nmax '+model_eval['title']+f' = {max_var:.2f}'
                print(stat_print_string)
                    # dist = 0.1*psc
                    # model_eval['ax'].text(0.05,0.4-dist,stat_print_string,transform=model_eval['ax'].transAxes,ha='left',va='center')
                    # psc+=1
    model_eval['ax'].set_xlabel(model_eval['title'])
    model_eval['ax'].set_ylim(ylim)
    model_eval['ax'].set_xlim(model_eval['xlim'])
    model_eval['ax'].set_xticks(model_eval['xticks'])
    model_eval['ax'].set_xticks(model_eval['xminorticks'],minor=True)
    model_eval['ax'].grid(alpha=1,color='lightgray',which='both')
    model_eval['ax'].text(0.05,0.95,atoz[count]+')',transform=model_eval['ax'].transAxes,va='top',ha='left')
    model_eval['ax'].axhline(y=0.001*(beach_list['TJRE']['r']-beach_list['PB']['r']),color='yellowgreen',linestyle='dashed',linewidth=1.0)
    if count>0:
        model_eval['ax'].set_ylabel('')
        model_eval['ax'].set_yticklabels([''])
    if count==0:
        model_eval['ax'].set_ylabel('y (km)')
    count+=1

# RMSE is only metric that wants to be logarithmic
if RMSE in model_eval_list:
    RMSE['ax'].set_xscale('log')

# this assumes WSS is your rightmost axis. Adjust accordingly if not true.
wqfun.add_beaches_to_ax_RHS(WSS['ax'],TJRE_line=True,xaxscale=0.05)

# need a legend
NRMSE['ax'].text(0.7,0.76,'1D',color=c1D,transform=NRMSE['ax'].transAxes,ha='left',bbox=props,fontsize=10)
NRMSE['ax'].text(0.7,0.69,'1DR',color=c1DR,transform=NRMSE['ax'].transAxes,ha='left',bbox=props,fontsize=10)
NRMSE['ax'].text(0.7,0.62,'1DC',color=c1DC,transform=NRMSE['ax'].transAxes,ha='left',bbox=props,fontsize=10)


# show or save plot
plt.subplots_adjust(bottom=0.2)
# plt.show(block=False)
# plt.pause(0.1)
home = '/Users/elizabethbrasseale/Projects/Water quality/'
outfn = home+'WQ_plots/paper figures/Figure_07_F.eps'
plt.savefig(outfn,format='eps')
