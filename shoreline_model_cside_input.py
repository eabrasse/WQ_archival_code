#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1D model of nearshore dye advection and loss
using alongshore varying velocity extracted from
COAWST model out to at 5-m isobath
and rotated to alongshore/crosshore coordinates
and interpolated to regular grid

"""

# setup
import os
import sys
alp = os.path.abspath('/Users/elizabethbrasseale/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun

import cmocean as cmo
import pickle
import numpy as np

from datetime import datetime, timedelta


home = '/Users/elizabethbrasseale/Projects/Water quality/'

###### SET DECAY SCALE PARAMETERS ######
#
# temporal decay is the rate at which bacteria die off
# this is the e drop scale
# 
# default is 10 days = 86400 secondds
time_decay_scale = 10.0 * 24 * 60 * 60

#
# Spatial decay scale is now hardcoded after tuning - EB 4/29/22
#
# spatial decay is used to determine a rate
# of offshore diffusion. Given a typical
# velocity, over what distance do you expect
# offshore diffusion to become signficant?
#
# Note: determined from 
#
# default is 10km = 10000 m
# length_decay_scale = 3750
# length_decay_scale = 7122.5

# if you want a shorter time series,
# define it here
t0 = 0
t1 = None

###### GET VELOCITY FROM WAVE DATA ######
#
cside_uv_fn = home + 'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Duv = pickle.load(open(cside_uv_fn,'rb'))
V = Duv['v_rot_int'][:]

############################################
###### ADVECTION DIFFUSION MODEL ######

###### SET UP TIME (t), SPACE (y) GRID ######
y = Duv['rshore'][:]

nt0,ny = Duv['dye_01'].shape

scale=200
nt = nt0*scale
t0 = np.linspace(0,nt0,nt0)
t = np.linspace(0,nt0,nt)

dt = 3600/scale
dy = np.diff(y)[0]

###### INTERPOLATE V ONTO NEW TIME GRID ######
#extend it over more time steps
Vmod_sm = np.zeros((nt,ny))
for tt in range(nt0):
    for j in range(ny):
        Vmod_sm[tt*scale:(tt+1)*scale,j] = V[tt,j]

    
###### SET UP INFLOW FORCING AT PUNTA BANDERA ######
y_ind = 89 # derived from CSIDE model; location of Punta Bandera on coastline
PB_in = 0.011 # EAB 4/29/22 result from PB_in tuning

###### DERIVE DECAY SCALES FROM USER INPUT PARAMETERS ######
kt = -1.0 / time_decay_scale # 1-day half-life, this number is negative ~ -8e-06
# kd = -1.0 * 0.1 / length_decay_scale # amount of diffusion for a typical velocity (10cm/s) and distance (10km)
kd = 1.3e-5

###### INITIALIZE DYE CONCENTRATION (c) ######
c = np.zeros((nt,ny))
adv = np.zeros((nt,ny-2))
c[0,y_ind] = PB_in

###### SOLVE ######
for tt in range(1,nt):

    adv_plus = Vmod_sm[tt-1,1:-1]*(c[tt-1,1:-1]-c[tt-1,0:-2])/dy
    adv_minus = Vmod_sm[tt-1,1:-1]*(c[tt-1,2:]-c[tt-1,1:-1])/dy
    adv[tt,:] = adv_plus
    adv[tt,Vmod_sm[tt-1,1:-1]<0] = adv_minus[Vmod_sm[tt-1,1:-1]<0]
    
    bac_decay = kt*c[tt-1,1:-1] #bacterial deacy should be negative because k is negative by c is positive
    offshore_diffusion = kd*c[tt-1,1:-1]
    
    c[tt,1:-1] = c[tt-1,1:-1] + dt*(-adv[tt,:] + bac_decay + offshore_diffusion)
    c[tt,y_ind] = PB_in

###### SUBSAMPLE c BEFORE SAVING ######
# necessary to save memory, python won't pickle a 4 GB file
c_ss = c[::scale,:] # interpolates it onto the model-derived, coarser 1-hr time grid
adv_ss = adv[::scale,:]
v_ss = Vmod_sm[::scale,:]


###### SAVE VARIABLES IN PICKLE FILE ######
D = {'c':c_ss,'t':t0,'y':y,'kt':kt,'kd':kd,'adv':adv_ss,'v':v_ss}
outfn = home+f'/WQ_data/adv_diff_model/CSIDE_recycled_input_tuned_subsampled_kd{kd:.2E}_PB_in{PB_in:0.3f}.p'
pickle.dump(D,open(outfn,'wb'))
