#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot 3-bin histograms comparing dye concentrations in 1D models with COAWST model
at 4 locations
The 1D model used here is 
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
import string
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
from matplotlib.ticker import LogFormatter

# load in some helpful plotting tools
c20 = plt.get_cmap('tab20b',20)
atoz = string.ascii_lowercase
    
plt.close('all')

# load in 1D model
home = '/Users/elizabethbrasseale/Projects/Water quality/'
model_name_list = ['CSIDE_SAsm100','U_isobath5m_sawc_autotune_kd']
D = wqfun.get_shoreline_models(model_name_list)
model0 = D['U_isobath5m_sawc_autotune_kd']
COAWST0 = D['CSIDE_SAsm100']

# get beach locations
beach_name_list = ['PB','PTJ','IB','SS','HdC']
beach_list = wqfun.get_beach_location(beach_name_list)
y = model0['y'][:]-beach_list['PB']['r']


# define cut-off for beach advisory condition
hazard = 5e-4

#generate figure and axis handles
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,fh))
gs = GridSpec(1,2)


# let's try something different
model = {}
COAWST = {}
model['dye'] = model0['dye'][:]+1e-12 #add 1e-12 because 0s are difficult in logspace
nt,ny = model0['dye'].shape
norm = 100/nt
axtickrange = range(0,101,20)
COAWST['dye'] = COAWST0['dye'][:]+1e-12 #add 1e-12 because 0s are difficult in logspace
count=0

dyemin = np.min([COAWST['dye'].min(),model['dye'].min()])
dyemax = np.max([COAWST['dye'].max(),model['dye'].max()])
logbin_edges = np.logspace(-4,-2,11) # this is the space we care about, but...
logbin_edges = np.insert(logbin_edges,0,1e-8)
logbin_edges = np.insert(logbin_edges,0,dyemin) # we want every value accounted for,
logbin_edges = np.append(logbin_edges,dyemax) # so include the data limits
logbins = 0.5*np.diff(logbin_edges) + logbin_edges[:-1]

yt = np.tile(np.reshape(y,(ny,1)),len(logbins))
lbt = np.transpose(np.tile(np.reshape(logbins,(len(logbins),1)),ny))

# colorscale plotting parameters
# dyemax = 5e-2
# dyemin = 1e-4
vmax = dyemax
vmin = dyemin
colmap = cmo.cm.matter
normal = matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax)

cmap = plt.get_cmap(cmo.cm.matter)
# cmap = cmo.tools.cmap('matter', N=len(logbins))

for data in model,COAWST:
    
    data['near_zero'] = np.sum(data['dye']<1e-8,axis=0)
    data['intermed'] = np.sum((data['dye']>1e-8)&(data['dye']<hazard),axis=0)
    data['big'] = np.sum(data['dye']>hazard,axis=0)
    
    data['hist'] = np.zeros((ny,len(logbins)))
    for yy in range(ny):
        hist,bin_edges = np.histogram(data['dye'][:,yy],bins=logbin_edges)
        data['hist'][yy,:] = hist*norm
    
    data['ax'] = fig.add_subplot(gs[count])
    
    for lbin in range(len(logbins)):
        left = np.sum(data['hist'][:,:lbin],axis=1)
        colind = lbin/len(logbins)
        data['ax'].barh(0.001*y,data['hist'][:,lbin],left=left,color=cmap(colind)) #this is not normed correctly, but it's a start)
    # p = data['ax'].pcolormesh(data['hist'][:,::-1],0.001*yt,lbt[:,::-1],cmap = colmap,norm=normal)

    data['ax'].plot(norm*data['near_zero'],0.001*y,color='k',linestyle='solid',lw=1.0)
    data['ax'].plot(norm*data['intermed']+norm*data['near_zero'],0.001*y,color='k',linestyle='solid',lw=1.0)
    data['ax'].plot(norm*data['big']+norm*data['intermed']+norm*data['near_zero'],0.001*y,color='k',linestyle='solid',lw=1.0)
    
    
    data['ax'].set_xticks(axtickrange)
    data['ax'].set_xticklabels([f'{r}'+r'%' for r in axtickrange])
    data['ax'].set_xlabel('Percent of time steps')
    data['ax'].set_ylabel('y (km)')
    if count>0:
        data['ax'].set_ylabel('')
        data['ax'].set_yticks([])
    data['ax'].set_ylim([0,0.001*y.max()])
    data['ax'].set_xlim([0,100])
    
    wqfun.add_beaches_to_ax_RHS(data['ax'],TJRE_line=True,xaxscale=0.03)
    
    count+=1
model['ax'].text(0.1,0.95,'a) 1D model',transform=model['ax'].transAxes,ha='left',va='top')
COAWST['ax'].text(0.1,0.95,'b) SD Bight model',transform=COAWST['ax'].transAxes,ha='left',va='top')

# add legend
model['ax'].text(0.48,0.1,r'0<$C$<1e-8',transform=model['ax'].transAxes,color='k',rotation=63)
model['ax'].text(0.66,0.05,r'1e-8<$C$<5e-4',transform=model['ax'].transAxes,color='k',rotation=70)
model['ax'].text(0.82,0.1,r'5e-4<$C$',transform=model['ax'].transAxes,color='k',rotation=75)

norm = matplotlib.colors.BoundaryNorm(logbin_edges[1:-1], cmap.N, extend='both')
cbaxes = inset_axes(COAWST['ax'], width="6%", height="60%", loc='center right',bbox_transform=COAWST['ax'].transAxes,bbox_to_anchor=(0.3,0.,1,1))
fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=cbaxes, orientation='vertical',ticks = [1e-8,1e-4,1e-3,1e-2],
             label="bins")#,format='%0.0e')
# override temporarily (probably a better way to do this)
# matplotlib.rcParams['mathtext.fontset'] = 'dejavusans'
# cbaxes.set_yticklabels([r'$10^{-8}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'],fontname='dejavusans')
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
# cbaxes.yaxis.set_major_formatter(LogFormatter())
cbaxes.set_yscale('log')

# show or save plot
# plt.tight_layout()
plt.subplots_adjust(bottom=0.15,right=0.8)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_11_E.eps'
plt.savefig(outfn,format='eps')
