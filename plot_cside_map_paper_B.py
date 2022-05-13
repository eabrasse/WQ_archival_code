#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot map of COAWST/CSIDE/SD Bight model with landmarks for context
"""

# setup
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
import wqfun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmo
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import string
from matplotlib.colors import LogNorm
s10 = plt.get_cmap('Set1')

plt.close('all')
fs_small=10
fs_big = 12


def add_features(ax,axes):
    # use cartopy features to plot landmarks on map
    ax.set_extent(axes, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    # ax.add_feature(cfeature.RIVERS)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'

    ax.add_feature(states_provinces, edgecolor='gray')

def plot_domain(lonr,latr,ax,col):
    # outline different model level domains
    # note: this program was more useful when I was plotting all 4 domains.
    # eventually, I only decided to plot the relevant LV4 domain
    lw = 1
    ls = 'solid'
    ax.plot(lonr[0,:],latr[0,:],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[-1,:],latr[-1,:],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[:,0],latr[:,0],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[:,-1],latr[:,-1],color=col,lw=lw,linestyle=ls)

props = dict(boxstyle='round', fc='white',ec='None',alpha=1.0)

# load in model domains
home = '/Users/elizabethbrasseale/Projects/Water quality/'
grid_dir = home+'WQ_data/'
# gridfn1 = grid_dir + 'GRID_SDTJRE_LV1_rx020.nc'
# dsg1 = nc.Dataset(gridfn1)
# gridfn2 = grid_dir + 'GRID_SDTJRE_LV2_rx020.nc'
# dsg2 = nc.Dataset(gridfn2)
# gridfn3 = grid_dir + 'GRID_SDTJRE_LV3_rx020.nc'
# dsg3 = nc.Dataset(gridfn3)
gridfn4 = grid_dir + 'GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
dsg4 = nc.Dataset(gridfn4)

# load in beach location information
beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
beach_list = wqfun.get_beach_location(beach_name_list)

# define axis ranges
# note: I previously included an inset axis showing all 4 domains. 
axes_big = [-125,-110,25,38]
axes_med = [-117.28,-117.05,32.4,32.78]
axes_small = [-117.23,-117.09,32.43,32.7]

# generate figure and name axis handles
fw,fh = wqfun.gen_plot_props(fs_big=fs_big,fs_small=fs_small)
fig = plt.figure(figsize=(2*fw,fh))
ax0 = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
ax1 = fig.add_subplot(1,2,2,projection=ccrs.PlateCarree())
axs = [ax0,ax1]

# add map features to the left axis
# add_features(ax=ax0,axes=axes_small)
ax0.set_extent(axes_small, crs=ccrs.PlateCarree())
ax1.set_extent(axes_small, crs=ccrs.PlateCarree())
bathy_lvls = [5,10,20]
lw = 0.4
# plot_domain(dsg4['lon_rho'][:],dsg4['lat_rho'][:],ax0,s10(4))

#load in location of wave buoy from wave buoy extraction dataset
wave_fn = grid_dir+'shoreline_variables_SZ_2017–2018.p'
Dbuoy = pickle.load(open(wave_fn,'rb'))

# ADD BAHTYMETRY
# bathy contours & wave buoy location in both axes
for ax in axs:
    ax.contour(dsg4['lon_rho'][:],dsg4['lat_rho'][:],dsg4['mask_rho'][:],colors='k',linewidths=lw)
    ax.contour(dsg4['lon_rho'][:],dsg4['lat_rho'][:],dsg4['h'][:],levels=bathy_lvls,colors='k',linewidths=lw)
    # plot wave buoy
    ax.plot(Dbuoy['buoylon'],Dbuoy['buoylat'],mfc='cornflowerblue',mec='k',marker='d',markersize=10)
ax0.text(Dbuoy['buoylon']-0.01,Dbuoy['buoylat']-0.028,'wave buoy',fontsize=fs_small,rotation=0,color='w',ha='center',va='center')

# colormesh + colorbar in only ax0
cmap = cmo.cm.deep
cmap.set_bad(color='w', alpha=1)
h = dsg4['h'][:]
h[dsg4['mask_rho'][:]<1]=np.nan
p0=ax0.pcolormesh(dsg4['lon_rho'][:],dsg4['lat_rho'][:],h,cmap=cmap,vmax=40,shading='Nearest')
cbaxes = inset_axes(ax0, width="6%", height="40%", loc='lower left',bbox_transform=ax0.transAxes,bbox_to_anchor=(-0.5,0,1,1))
cb = fig.colorbar(p0, cax=cbaxes, orientation='vertical',ticks=[0,5,10,20,40])
cb.set_label('Seafloor depth (m)',fontsize=fs_small)
cbaxes.invert_yaxis()

# highlight shoreline
PBind = np.argmin(np.abs(Dbuoy['rshore']-beach_list['PB']['r']))
ax0.scatter(Dbuoy['lonshore'][PBind:],Dbuoy['latshore'][PBind:],c='red',s=5,marker=0)
for lvl in bathy_lvls:
    cbaxes.plot([0,1],[lvl]*2,'k',linewidth=lw)

# add an inset axis
ax_in = inset_axes(ax0,height='100%',width='100%',loc='upper right',\
bbox_to_anchor=(-0.7, .6, .7, .7),bbox_transform=ax0.transAxes,\
axes_class=cartopy.mpl.geoaxes.GeoAxes,axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
add_features(ax=ax_in,axes=axes_big)
# plot_domain(dsg4['lon_rho'][:],dsg4['lat_rho'][:],ax_in,s10(4))
ax_in.plot(dsg4['lon_rho'][:].mean(),dsg4['lat_rho'][:].mean(),mfc=s10(4),mec='k',marker='*',markersize=10)

#ADD SECOND PLOT OF DYE
surface_dye_fn = home+'WQ_data/extractions2017/surface_dye_01_July.p'
D = pickle.load(open(surface_dye_fn,'rb'))
p1=ax1.pcolormesh(dsg4['lon_rho'][:],dsg4['lat_rho'][:],D['Jul11']['dye_01'][:],norm=LogNorm(vmin=1e-4,vmax=5e-2),cmap=cmo.cm.matter,shading='nearest')
cbaxes = inset_axes(ax1, width="6%", height="40%", loc='lower left',bbox_transform=ax1.transAxes,bbox_to_anchor=(0.,0,1,1))
cb = fig.colorbar(p1, cax=cbaxes, orientation='vertical')

# ADD BEACHES TO BOTH PLOTS
beach_list['TJRE']['lon'] = beach_list['TJRE']['lon']+0.02
beach_list['PB']['lon'] = beach_list['PB']['lon']+0.005
for beach in beach_name_list:
    if beach=='PB':
        beach_color='m'
        mark = '<'
    elif beach=='TJRE':
        beach_color='yellowgreen'
        mark = '<'
    else:
        beach_color='yellow'
        mark = 'o'
    for ax in ax0,ax1:
        ax.plot(beach_list[beach]['lon'],beach_list[beach]['lat'],marker=mark,mec='k',mfc=beach_color,markersize=8)
    
    ax0.text(beach_list[beach]['lon']+0.01,beach_list[beach]['lat'],beach,va='center',ha='left',color='k',bbox=props,fontsize=fs_small)

#add labels
ax0.text(0.95,0.98,'b)',ha='right',va='top',transform=ax0.transAxes)
ax0.set_xticks([-117.2,-117.1])
ax0.set_xlabel('Longitude')
# ax0.set_yticks([32.45,32.6,32.75])
ax0.set_ylabel('Latitude')
ax0.yaxis.set_label_position("right")
ax0.yaxis.tick_right()
ax0.set_yticks([32.5,32.6])

ax1.text(0.95,0.98,'c)',ha='right',va='top',transform=ax1.transAxes)
ax1.get_xaxis().set_visible(True)
ax1.set_xlabel('Longitude')
ax1.set_xticks([-117.2,-117.1])
ax1.get_yaxis().set_visible(True)
ax1.set_ylabel('Latitude')
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.set_yticks([32.5,32.6])

ax_in.text(0.05,0.05,'a)',ha='left',va='bottom',transform=ax_in.transAxes)
ax_in.set_xlabel('Longitude',fontsize=10)
ax_in.set_ylabel('Latitude',fontsize=10)
ax_in.set_xticks([-120,-115])
ax_in.set_yticks([30,35])
ax_in.tick_params(length=3)

# show or save plot
plt.subplots_adjust(bottom=0.13,left=0.2,wspace=0)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'WQ_plots/paper figures/Figure_01_C.jpg'
plt.savefig(outfn,format='jpg',dpi=600)