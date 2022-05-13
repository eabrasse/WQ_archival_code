#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process velocity (rotate and interpolate) output of CSIDE extraction
"""

# setup
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
alp = os.path.abspath('/Users/elizabethbrasseale/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmo
c10 = plt.get_cmap('tab10')

plt.close('all')

def smooth(vec,window_width):
    half_window = int(window_width/2)
    # cumsum_vec = np.cumsum(np.insert(vec, 0, 0))
    # vec_padded = np.pad(vec, (half_window, window_width-1-half_window), mode='edge')
    vec_padded = np.pad(vec, (half_window, window_width-1-half_window), mode='constant',constant_values=(np.mean(vec[:10]),np.mean(vec[-10:])))
    cumsum_vec = np.cumsum(np.insert(vec_padded,0,0))
    new_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return new_vec

home = '/Users/elizabethbrasseale/Projects/Water Quality/'
uv_fn = home + 'WQ_data/extractions2017/shoreline_uv_05m.p'
D=pickle.load(open(uv_fn,'rb'))

# var_list = ['dye_01','dye_02','ot','lon_rho','lat_rho','mask_rho','shorelon','shorelat']
# var_list = ['dye_01','dye_02','ot','lon_rho','lat_rho','mask_rho','iis','jjs','iib','jjb','Dwave','Hwave','Lwave','zeta']

for var in D.keys():
    locals()[var]=D[var]

nt,nj = u.shape

lonshore = np.zeros((nj))
latshore = np.zeros((nj))
for j in range(nj):
    lonshore[j] = lon_rho[jjs[j],iis[j]]
    latshore[j] = lat_rho[jjs[j],iis[j]]

dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

x_rho,y_rho = zfun.ll2xy(lon_rho,lat_rho,lon_rho.min(),lat_rho.min())
xshore,yshore = zfun.ll2xy(lonshore,latshore,lon_rho.min(),lat_rho.min())
# x10m,y10m = zfun.ll2xy(lon10m,lat10m,lon_rho.min(),lat_rho.min())
dxs = np.diff(xshore)
dys = np.diff(yshore)
drs = np.sqrt(dxs**2 + dys**2)
rshore = np.cumsum(drs)
rshore = np.insert(rshore,0,0,axis=0)
nr = nj


#new grid is evenly spaced
rshore_new = np.linspace(0,rshore[-1],nr)

#calculate along- and cross-shore velocity
grid_angle = np.mean(np.arctan2(np.diff(lat_rho,axis=1),np.diff(lon_rho,axis=1)))
# V = np.sqrt(u**2+v**2)
w0 = u + 1j*v
w = w0 * np.exp(1j*grid_angle)

theta_SA = (shoreangle)*np.pi/180
# uu_SA = u*np.cos(theta_SA)+v*np.sin(theta_SA)
# vv_SA = -u*np.sin(theta_SA)+v*np.cos(theta_SA)
w_SA = w * np.exp(1j*theta_SA)
u_SA = np.real(w_SA)
v_SA = np.imag(w_SA)

window_width=50
# half_window = int(window_width/2)
# cumsum_vec = np.cumsum(np.insert(theta_PA, 0, 0))
# cumsum_vec_padded = np.pad(cumsum_vec, (half_window, window_width-1-half_window), mode='edge')
# theta_PA = (cumsum_vec_padded[window_width:] - cumsum_vec_padded[:-window_width]) / window_width
theta_SA = smooth(theta_SA,window_width)

# uu_PA = u*np.cos(theta_PA)+v*np.sin(theta_PA)
# vv_PA = -u*np.sin(theta_PA)+v*np.cos(theta_PA)
uu_new = np.zeros((nt,nr))
vv_new = np.zeros((nt,nr))
for t in range(nt):
    uu_new[t,:] = np.interp(rshore_new,rshore,u_SA[t,:])
    vv_new[t,:] = np.interp(rshore_new,rshore,v_SA[t,:])

dye_fn = home + 'WQ_data/extractions2017/shoreline_dye_waves_05m_interp.p'
Ddye = pickle.load(open(dye_fn,'rb'))

Dout = {}
var_list = ['dye_01','dye_02','Dwave','Lwave','Hwave','zeta']
for var in var_list:
    Dout[var] = Ddye[var]

Dout['rshore'] = rshore_new
Dout['rshore_old'] = rshore
Dout['u'] = u
Dout['v'] = v
Dout['u_rot'] = u_SA
Dout['v_rot'] = v_SA
Dout['u_rot_int'] = uu_new
Dout['v_rot_int'] = vv_new
Dout['ot'] = ot
Dout['shoreangle'] = shoreangle

outfn = home + 'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
pickle.dump(Dout,open(outfn,'wb'))

# fig = plt.figure(figsize=(12,10))
# gs = GridSpec(2,2)
# axu_SA = fig.add_subplot(gs[0,0])
# axu_saved = fig.add_subplot(gs[0,1])
# axv_SA = fig.add_subplot(gs[1,0])
# axv_saved = fig.add_subplot(gs[1,1])

# saved_fn = home + 'WQ_data/extractions2017/shoreline_uv_05m_interp_sm100_SA.p'
# Dsaved=pickle.load(open(saved_fn,'rb'))
# u_saved = Dsaved['u'][:]
# v_saved = Dsaved['v'][:]
#
# umax = 0.75*np.max(np.abs([u_SA.min(),u_SA.max()]))
# vmax = 0.75*np.max(np.abs([v_SA.min(),v_SA.max()]))
# p = axv_SA.pcolormesh(dt_list,0.001*rshore,np.transpose(v_saved-v_SA),vmin=-vmax,vmax=vmax,cmap='RdBu_r',shading='nearest')
# axv_SA.set_ylabel('Alongshore distance (km)',fontweight='bold')
# axv_SA.set_title('saved v -rotated v')
# cbaxes = inset_axes(axv_SA, width="4%", height="80%", loc=4,bbox_transform=axv_SA.transAxes,bbox_to_anchor=(-0.075,0.,1,1))
# cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
# cb.set_label('velocity')
#
# p = axu_SA.pcolormesh(dt_list,0.001*rshore,np.transpose(u_saved-u_SA),vmin=-vmax,vmax=vmax,cmap='RdBu_r',shading='nearest')
# axu_SA.set_ylabel('Alongshore distance (km)',fontweight='bold')
# axu_SA.set_title('saved u - rotated u')
#
# p = axv_saved.pcolormesh(dt_list,0.001*rshore_new,np.transpose(v_saved-vv_new),vmin=-vmax,vmax=vmax,cmap='RdBu_r',shading='nearest')
# axv_saved.set_ylabel('Alongshore distance (km)',fontweight='bold')
# axv_saved.set_title('saved v - rotated and interpolated v')
#
# p = axu_saved.pcolormesh(dt_list,0.001*rshore_new,np.transpose(u_saved-uu_new),vmin=-vmax,vmax=vmax,cmap='RdBu_r',shading='nearest')
# axu_saved.set_ylabel('Alongshore distance (km)',fontweight='bold')
# axu_saved.set_title('saved u - rotated and interpolated u')
#
# for ax in axu_SA,axu_saved:
#     ax.set_xticklabels([''])
#     ax.get_xaxis().set_visible('False')
# for ax in axv_SA,axv_saved:
#     ax.set_xlabel('Time',fontweight='bold')
#     ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
#     plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
#
#
# plt.tight_layout()
# plt.show(block=False)
# plt.pause(0.1)
