#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

home = '/data0/ebrasseale/'

dir0 = '/data0/NADB2017/NADB2017_0_NEW/'

fn0 = dir0+'ocean_his_NADB_0_new_00021.nc'
fn1 = dir0+'ocean_his_NADB_0_new_00022.nc'

ds0 = nc.Dataset(fn0)
ds1 = nc.Dataset(fn1)

ot0 = ds0['ocean_time'][:]
td00 = datetime(2017,7,8,12)-datetime(1999,1,1)
td01 = datetime(2017,7,9,12)-datetime(1999,1,1)
ind00 = np.argmin(np.abs(ot0-td00.total_seconds()))
ind01 = np.argmin(np.abs(ot0-td01.total_seconds()))

ot1 = ds1['ocean_time'][:]
td10 = datetime(2017,7,10,12)-datetime(1999,1,1)
td11 = datetime(2017,7,11,12)-datetime(1999,1,1)
ind10 = np.argmin(np.abs(ot1-td10.total_seconds()))
ind11 = np.argmin(np.abs(ot1-td11.total_seconds()))

dye_01_0 = ds0['dye_01'][:]
dye_01_1 = ds1['dye_01'][:]

Jul8 = {}
Jul8['dye_01'] = ds0['dye_01'][ind00,-1,:,:]
Jul8['ocean_time'] = ot0[ind00]

Jul9 = {}
Jul9['dye_01'] = ds0['dye_01'][ind01,-1,:,:]
Jul9['ocean_time'] = ot0[ind01]

Jul10 = {}
Jul10['dye_01'] = ds1['dye_01'][ind10,-1,:,:]
Jul10['ocean_time'] = ot1[ind10]

Jul11 = {}
Jul11['dye_01'] = ds1['dye_01'][ind11,-1,:,:]
Jul11['ocean_time'] = ot1[ind11]


D = dict()
D['Jul8'] = Jul8
D['Jul9'] = Jul9
D['Jul10'] = Jul10
D['Jul11'] = Jul11

outfn = home + 'WQ_data/surface_dye_01_July.p'
pickle.dump(D,open(outfn,'wb'))

ds0.close()
ds1.close()