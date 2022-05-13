#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code plots velocity from the COAWST model and the 1D models
as a time series and as a scatter plot to demonstrate agreement
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

props = dict(boxstyle='round', facecolor='white',edgecolor='white')
# plt.rcParams['text.usetex'] = True
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
# plt.rcParams["font.family"] = ["Latin Modern Roman"]
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
    
# c10 = plt.get_cmap('tab10',10)
# generate colors from colormap
c20c = plt.get_cmap('tab20c',20)
# cC = 'k'
c1D = c20c(0)
# c1DC = c20c(6)
c1DR = c20c(9)
c_list = [c1D,c1DR]


# load in files
plt.close('all')
home = '/Users/elizabethbrasseale/Projects/Water quality/'

#I use "cside" and "coawst" interchangeably - they both refer to the SD Bight model
# CSIDE was the name of a larger combined field-modeling effort to understand
# dynamics in this region. The lab group I joined used both phrases interchangeably
# so I went back and forth on notation. Just to be confusing, it's referred to in the
# manuscript as the SD Bight model.

# velocities from the COAWST model: extracted, rotated, and interpolated
cside_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Dcside = pickle.load(open(cside_fn,'rb'))
v_COAWST = np.mean(Dcside['v_rot_int'][1:],axis=1)
# build up a list of datetimes to be used as the x-axis for timeseries plots
ot=Dcside['ot'][1:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))
nt = ot.shape[0]

# load in wave variables from wave buoy to calculate velocities
wave_fn = home+'WQ_data/shoreline_variables_SZ_2017–2018.p'
Dbuoy = pickle.load(open(wave_fn,'rb'))
# we're using velocities at the 5m isobath, but the code I wrote calculates both.
# So we're loading in both, but only using one
WaveBuoy,Isobath5m = wqfun.wavebuoy_to_5m_isobath(Dbuoy)

v1D = Isobath5m['V'][:nt]
v1Dr = Isobath5m['Sxy'][:nt]/Isobath5m['linear_fit']

# generate figure and axis handles
fw,fh = wqfun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,2*fh))
gs = GridSpec(10,2)
ax_vel = fig.add_subplot(gs[:2,:])
ax = fig.add_subplot(gs[2:4,:])
ax1 = fig.add_subplot(gs[4:6,:])
ax0 = fig.add_subplot(gs[7:,0])
ax2 = fig.add_subplot(gs[7:,1])

#set linewidth
lw = 1.2

# load in beach properties. We only really need PB, because
# y-axis is measured in km from PB
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)

# plot alongshore-varying alongshore velocity from COAWST model
y = Dcside['rshore'][1:-1]-beach_list['PB']['r']
vmax = 0.75*np.max(np.abs([Dcside['v_rot_int'].min(),Dcside['v_rot_int'].max()]))
vmin = -vmax
velcol = 'RdBu_r'
velnorm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
bb2a = (0.1,0.35,1,1)
pv = ax_vel.pcolormesh(dt_list,0.001*y,np.transpose(Dcside['v_rot_int'][1:,1:-1]),norm=velnorm,cmap=velcol,shading='nearest')
ax_vel.text(0.05,0.95,'a)',transform=ax_vel.transAxes,ha='left',va='top')
ax_vel.set_ylabel('y (km)')
cbaxes = inset_axes(ax_vel, width="30%", height="8%", loc=6,bbox_transform=ax_vel.transAxes,bbox_to_anchor=bb2a)
cb = fig.colorbar(pv, cax=cbaxes, orientation='horizontal')
cbaxes.set_xlabel('m/s')
ax_vel.set_xticklabels([''])
ax_vel.set_xlabel('')

wqfun.add_beaches_to_ax_RHS(ax_vel,TJRE_line=False,xaxscale=0.02)

# now plot line plots of alongshore-mean alongshore velocity from COAWST and estimated (axA)
# and scatter plots (axB)
# match x axis (time axis) limits with the first plot 
xlim = ax_vel.get_xlim()
add_linreg = True # adds a linear regression
velname = ['v1D','v1DR']
x = np.linspace(-1,1,100)
count =0
for [axA,axB,v] in [[ax,ax0,v1D],[ax1,ax2,v1Dr]]:
    axA.plot(dt_list,v,color=c_list[count],lw=lw,label='calc from waves')
    axA.plot(dt_list,v_COAWST,color='k',lw=lw*0.5,linestyle='solid',label='COAWST')
    axA.set_ylabel(r'Velocity (m/s)')
    axA.grid(color='lightgray',alpha=1.0)
    axA.set_ylim([-0.75,0.75])
    axA.set_xlim(xlim)
    # axA.set_aspect(100)
    
    axB.axis([-0.5,0.5,-0.75,0.75])
    axB.set_xticks([-0.5,0,0.5])
    axB.set_xticks([-0.5,-0.25,0,0.25,0.5],minor=True)
    axB.set_aspect(1)
    axB.grid(color='lightgray',alpha=1.0,which='both')
    axB.scatter(v_COAWST,v,color=c_list[count],s=0.5)
    axB.set_xlabel(r'$\bar{v}_{C}$ (m/s)')
    axB.plot(x,x,color='m',linestyle='solid',lw=lw)
    
    if add_linreg:
        p = np.polyfit(v_COAWST,v,1)
        y = p[0]*x+p[1]
        axB.plot(x,y,linestyle='dashed',lw=2*lw,color='k')
        stat_string = velname[count]+' stats\n'+f'slope = {p[0]:}'+'\n'+f'intercept = {p[1]:}'+'\n\n'
        print(stat_string)
    
    count+=1


# Manually add a legend
ax.text(0.275,0.81,r'$\bar{v}_{C}$',color='k',transform=ax.transAxes,va='top',ha='left',bbox=props)
ax.text(0.275,0.95,r'$v_{1\mathrm{D}}$',color=c1D,transform=ax.transAxes,va='top',ha='left',bbox=props)
ax1.text(0.275,0.81,r'$\bar{v}_{C}$',color='k',transform=ax1.transAxes,va='top',ha='left',bbox=props)
ax1.text(0.275,0.95,r'$v_{1\mathrm{DR}}$',color=c1DR,transform=ax1.transAxes,va='top',ha='left',bbox=props)

# add axis labels
ax.set_xticklabels([''])
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d"))
plt.setp( ax1.xaxis.get_majorticklabels(), rotation=15, ha="right",rotation_mode='anchor')
ax1.set_xlabel('Date of 2017',fontweight='normal')
ax.text(0.02,0.98,'b)',transform=ax.transAxes,ha='left',va='top')
ax1.text(0.02,0.98,'c)',transform=ax1.transAxes,ha='left',va='top')
ax0.text(0.02,0.98,'d)',transform=ax0.transAxes,ha='left',va='top')
ax2.text(0.02,0.98,'e)',transform=ax2.transAxes,ha='left',va='top')

ax0.set_ylabel(r'$v_{1\mathrm{D}}$ (m/s)')
ax2.set_ylabel(r'$v_{1\mathrm{DR}}$ (m/s)')

# show or save plot
plt.subplots_adjust(top=0.98,bottom=0.08,hspace=0.08)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_03_G.jpg'
plt.savefig(outfn,format='jpg',dpi=600)
