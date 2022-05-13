#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1D model of nearshore dye advection and loss
using uniform velocity calcualted from
wave properties determined at 5-m isobath
from waves at an offshore location

"""

# setup
import os
import sys
alp = os.path.abspath('/Users/elizabethbrasseale/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import wqfun
import cmocean as cmo
import pickle
import numpy as np
import time
from datetime import datetime, timedelta

tic = time.perf_counter()
home = '/Users/elizabethbrasseale/Projects/Water quality/'

###### SET DECAY SCALE PARAMETERS ######
#
# temporal decay is the rate at which bacteria die off
# this is the e drop scale
# 
# default is 10 days = 86400 seconds - decay for norovirus
time_decay_scale = 10 * 24 * 60 * 60

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
t1 = 8977

###### GET VELOCITY FROM WAVE DATA ######
Rayleigh=False
#
# use wave data from "x_cside_vars_surfzone_calc.py"
wave_data = home+'WQ_data/shoreline_variables_SZ_2017–2018.p'

D = pickle.load(open(wave_data,'rb'))

WaveBuoy,Isobath5m = wqfun.wavebuoy_to_5m_isobath(D)
if Rayleigh:
    V = Isobath5m['Sxy']/Isobath5m['linear_fit']
else:
    V = Isobath5m['V'][:]


############################################
###### ADVECTION DIFFUSION MODEL ######

###### SET UP TIME (t), SPACE (y) GRID ######
# ultimately, we want to compare with exracted/processed COAWST data
# so start with those dimensions
shoreline_dye_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_05m_interp.p'
Ddye = pickle.load(open(shoreline_dye_fn,'rb'))
nt0,ny = Ddye['dye_01'].shape

# build time axis
scale=200
nt = nt0*scale
t0 = np.linspace(0,nt0,nt0)
t = np.linspace(0,nt0,nt)
y = np.linspace(0,Ddye['rshore'][-1],ny)

dt = 3600/scale
dy = np.diff(y)[0]

###### INTERPOLATE V ONTO NEW TIME GRID ######
Vmod = np.zeros((nt))

for tt in range(nt0):
    Vmod[tt*scale:(tt+1)*scale] = V[tt]

###### SET UP INFLOW FORCING AT PUNTA BANDERA ######
# these values are derived from the model grid tranformed
# from lat/lon to m in distance PB
y_ind = 89
if Rayleigh:
    PB_in = 0.01 # value derived from tuning in "shoreline_model_tune_PB_in.py"
else:
    PB_in = 0.008 # value derived from tuning in "shoreline_model_tune_PB_in.py"

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
    if Vmod[tt-1]>=0:
        adv = Vmod[tt-1]*(c[tt-1,1:-1]-c[tt-1,0:-2])/dy
    elif Vmod[tt-1]<0:
        adv = Vmod[tt-1]*(c[tt-1,2:]-c[tt-1,1:-1])/dy    
    bac_decay = kt*c[tt-1,1:-1] #bacterial deacy should be negative because k is negative by c is positive
    offshore_diffusion = kd*c[tt-1,1:-1]
    
    c[tt,1:-1] = c[tt-1,1:-1] + dt*(-adv + bac_decay + offshore_diffusion)
    c[tt,y_ind] = PB_in

###### SUBSAMPLE c BEFORE SAVING ######
# necessary to save memory, python won't pickle a 4 GB file
c_ss = c[::scale,:] # subsample onto the coarser model-derived, 1-hr time grid
v_ss = Vmod[::scale]

###### SAVE VARIABLES IN PICKLE FILE ######
D = {'c':c_ss,'t':t0,'y':y,'kt':kt,'kd':kd,'v':v_ss}
if Rayleigh:
    outfn = home+f'/WQ_data/adv_diff_model/CSIDE_uniformV_5miso_rayleigh_kd{kd:.2E}_tuned_shorenormal_PB_in{PB_in:0.3f}.p'
else:
    outfn = home+f'/WQ_data/adv_diff_model/CSIDE_uniformV_5miso_kd{kd:.2E}_tuned_shorenormal_PB_in{PB_in:0.3f}.p'
pickle.dump(D,open(outfn,'wb'))

toc = time.perf_counter()
model_time = f"shoreline model took {toc-tic:0.4f} seconds"
print(model_time)
