#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import wqfun
import cmocean as cmo
import pickle
import numpy as np
import time
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt

plt.close('all')

def solve_AV(nt,ny,Vmod,PB_in,kt,kd):
    ###### INITIALIZE DYE CONCENTRATION (c) ######
    c = np.zeros((nt,ny))
    adv = np.zeros((nt,ny-2))
    # c[0,:] = PB_in
    c[0,y_ind] = PB_in

    ###### SOLVE ######
    for tt in range(1,nt):
        adv_plus = Vmod[tt-1,1:-1]*(c[tt-1,1:-1]-c[tt-1,0:-2])/dy
        adv_minus = Vmod[tt-1,1:-1]*(c[tt-1,2:]-c[tt-1,1:-1])/dy
        adv[tt,:] = adv_plus
        adv[tt,Vmod[tt-1,1:-1]<0] = adv_minus[Vmod[tt-1,1:-1]<0]
        
        bac_decay = kt*c[tt-1,1:-1] #bacterial deacy should be negative because k is negative by c is positive
        offshore_diffusion = -kd*c[tt-1,1:-1]

        c[tt,1:-1] = c[tt-1,1:-1] + dt*(-adv[tt,:] + bac_decay + offshore_diffusion)
        c[tt,y_ind] = PB_in
    # necessary to save memory, python won't pickle a 4 GB file
    c_ss = c[::scale,:] # interpolates it onto the model-derived, coarser 1-hr time grid
    v_ss = Vmod[::scale]

    return c_ss,v_ss

def plot_PB_in_vs_wss(PB_in_array,wss_array):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel('PB_in')
    ax.set_ylabel('wss')
    ax.scatter(PB_in_array,wss_array,c='k',s=5)
    #find best two guesses and split the difference to make the next guess
    best_guess_ind = np.argmax(wss_array)
    wss_array_next_best = wss_array[np.arange(len(wss_array))!=best_guess_ind]
    PB_in_array_next_best = PB_in_array[np.arange(len(PB_in_array))!=best_guess_ind]
    next_best_guess_ind = np.argmax(wss_array_next_best)
    
    if PB_in_array[best_guess_ind]>PB_in_array_next_best[next_best_guess_ind]:
        ha_bg = 'left'
        ha_nbg = 'right'
    else:
        ha_bg = 'right'
        ha_nbg = 'left'

    ax.plot(PB_in_array[best_guess_ind],wss_array[best_guess_ind],mec='k',mfc='yellow',marker='*',markersize=10)
    ax.text(PB_in_array[best_guess_ind],wss_array[best_guess_ind]*0.9,'best guess\n'+f'PB_in={PB_in_array[best_guess_ind]:.4f}',color='k',va='top',ha=ha_bg)

    ax.plot(PB_in_array_next_best[next_best_guess_ind],wss_array_next_best[next_best_guess_ind],mec='k',mfc='gray',marker='*',markersize=10)
    ax.text(PB_in_array_next_best[next_best_guess_ind],wss_array_next_best[next_best_guess_ind]*0.9,'next best guess\n'+f'PB={PB_in_array_next_best[next_best_guess_ind]:.4f}',color='k',va='top',ha=ha_nbg)

    WSS_perc_diff = (wss_array[best_guess_ind]-wss_array_next_best[next_best_guess_ind])/wss_array_next_best[next_best_guess_ind]

    next_guess_diff = 0.5*np.abs(PB_in_array[best_guess_ind]-PB_in_array_next_best[next_best_guess_ind])
    next_guess = np.min([PB_in_array[best_guess_ind],PB_in_array_next_best[next_best_guess_ind]])+next_guess_diff

    ax.axvline(x=next_guess,color='k',linestyle='dashed')

    # ax.set_xscale('log')

    plt.show(block=False)
    plt.pause(0.1)

    return next_guess, WSS_perc_diff
#need location of PB
# beach_name_list = ['PB']
# beach_list = wqfun.get_beach_location(beach_name_list)
# PB_loc = beach_list['PB']['r']

# I want to run the model to automatically tune kd and then C0

home = '/Users/elizabethbrasseale/Projects/Water quality/'

###### SET DECAY SCALE PARAMETERS ######
# default is 10 days = 86400 seconds
time_decay_scale = 10 * 24 * 60 * 60
kt = -1.0 / time_decay_scale # 1-day half-life, this number is negative ~ -8e-06


###### GET VELOCITY FROM WAVE DATA ######

wave_data = home+'WQ_data/shoreline_variables_SZ_2017–2018.p'

D = pickle.load(open(wave_data,'rb'))


WaveBuoy,Isobath5m = wqfun.wavebuoy_to_5m_isobath(D)
#use kd from SD Bight model
kd = 1.3e-5

cside_uv_fn = home + 'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
Duv = pickle.load(open(cside_uv_fn,'rb'))
V = Duv['v_rot_int'][:]

PB_in0 = 0.01    

PB_in_array = np.array([PB_in0, PB_in0+0.005, PB_in0-0.005])

###### SET UP INFLOW FORCING AT PUNTA BANDERA ######
# PB_in = np.zeros((ny))
y_ind = 89 # derived from CSIDE model; location of Punta Bandera on coastline

# PB_in = 0.03 # average concentration at mouth in CSIDE model


# if you want a shorter time series,
# define it here
t0 = 0
t1 = 8977
# t1 = 50


############################################
###### ADVECTION DIFFUSION MODEL ######

###### SET UP TIME (t), SPACE (y) GRID ######
# nt0, ny = dye_01[t0:t1,:].shape
# shoreline_dye_fn = home+'WQ_data/extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
# Ddye = pickle.load(open(shoreline_dye_fn,'rb'))
nt0,ny = Duv['dye_01'][:].shape

scale=200
nt = nt0*scale
t0 = np.linspace(0,nt0,nt0)
t = np.linspace(0,nt0,nt)
# warning: hard coding ahead
# these values are derived from the model grid tranformed
# from lat/lon to m in distance from TJRE
y = np.linspace(0,Duv['rshore'][-1],ny)

dt = 3600/scale
dy = np.diff(y)[0]

###### INTERPOLATE V ONTO NEW TIME GRID ######
Vmod = np.zeros((nt,ny))
for tt in range(nt0):
    for j in range(ny):
        Vmod[tt*scale:(tt+1)*scale,j] = V[tt,j]

wss_array = np.array([])
#using the same V, solve using different kd
for PB_in in PB_in_array:
    print(f'Solving shoreline model with PB_in = {PB_in:0.4f}...')
    tic = time.perf_counter()
    c,v=solve_AV(nt,ny,Vmod,PB_in,kt,kd)
    toc = time.perf_counter()
    model_time = f"shoreline model took {toc-tic:0.4f} seconds"
    print(model_time)
    # evaluate wss by shoreline location, so that you're comparing with the mean
    # dye at that location and not the model mean
    wss_alongshore = np.zeros((ny-y_ind))
    for j in range(ny-y_ind):
        wss_alongshore[j] = wqfun.willmott(c[:,y_ind+j],Duv['dye_01'][:,y_ind+j])
    wss = np.mean(wss_alongshore)
    wss_array= np.append(wss_array,[wss],axis=0)


next_guess, WSS_perc_diff=plot_PB_in_vs_wss(PB_in_array,wss_array)
print(f'WSS improvement fraction {WSS_perc_diff:.8f}')

while True:
    while True:
        yn = input(f'Try next guess PB_in = {next_guess:0.4f}(dashed line)?\n(enter y for yes, o for override value, n for exit & save, or x for exit without saving)\n')
        if yn in ['y','n','o','x']:
            break
        else:
            continue
    if yn in ['n','x']:
        break
    if yn=='y':
        PB_in_guess = next_guess
    elif yn=='o':
        while True:
            my_guess_str = input('Enter override PB_in guess (positive value, supports scientific notation, e.g. 2.5e-7)\n')
            try:
                my_guess = float(my_guess_str)
            except ValueError:
                print('You must enter a numeric value\n')
                continue
            if my_guess>0:
                PB_in_guess = my_guess
                break
            else:
                print('Value of PB_in must be positive\n')
                continue

    PB_in_array= np.append(PB_in_array,[PB_in_guess],axis=0)

    print(f'Solving shoreline model with PB_in = {PB_in_guess:0.4E}...')
    tic = time.perf_counter()
    c,v=solve_AV(nt,ny,Vmod,PB_in_guess,kt,kd)
    toc = time.perf_counter()
    model_time = f"shoreline model took {toc-tic:0.4f} seconds"
    print(model_time)
    
    # wss = wqfun.willmott(c[:,y_ind:],Ddye['dye_01'][:,y_ind:])
    wss_alongshore = np.zeros((ny-y_ind))
    for j in range(ny-y_ind):
        wss_alongshore[j] = wqfun.willmott(c[:,y_ind+j],Duv['dye_01'][:,y_ind+j])
    wss = np.mean(wss_alongshore)
    wss_array= np.append(wss_array,[wss],axis=0)
    
    plt.close('all')
    next_guess, WSS_perc_diff=plot_PB_in_vs_wss(PB_in_array,wss_array)
    print(f'WSS improvement fraction {WSS_perc_diff:.8f}')

if yn=='x':
    print('Exiting without saving')
else:
    best_guess_ind = np.argmax(wss_array)
    PB_in = PB_in_array[best_guess_ind]
    ###### SAVE VARIABLES IN PICKLE FILE ######
    D = {'c':c,'t':t0,'y':y,'kt':kt,'kd':kd,'v':v,'PB_in':PB_in}
    model_type_string = 'recycled'
    model_name = 'autotuned_'+model_type_string+f'_kd{kd:.4E}_AVv_5miso_PB_in{PB_in:0.3f}'
    data_outfn = home+f'WQ_data/adv_diff_model/'+model_name+'.p'
    print(f'Saving to {data_outfn}')
    pickle.dump(D,open(data_outfn,'wb'))
    fig_outfn = home+f'WQ_plots/adv diff model/'+model_name+'.png'
    plt.savefig(fig_outfn)
    plt.close()

