#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experime.
"""

# setup

import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc


def get_shore_inds(grid_fn):
    
    ds = nc.Dataset(grid_fn)
    var_list = ['mask_rho','lon_rho','lat_rho']
    # I think locals doesn't work in a module
    # for var in var_list:
    #     locals()[var] = ds[var][:]
    mask_rho = ds['mask_rho'][:]
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]

    ny,nx = lat_rho.shape
    iis = np.zeros((ny),dtype=int)
    jjs = range(ny)
    lonshore = np.zeros((ny))
    latshore = np.zeros((ny))

    #extract shoreline location
    for j in jjs:
    
        # find the edge of the mask
        mask_diff = np.where(np.diff(mask_rho[j,:]))[0]
    
        #if multiple edges, north of TJRE
        if (len(mask_diff)>1)&(lat_rho[j,0]>32.6):
            #look for the edge closest to the previously identified edge
            iis[j] = mask_diff[np.argmin(np.abs(iis[j-1]-mask_diff))]
        
        #if multiple edges, south of TJRE
        elif (len(mask_diff)>1)&(lat_rho[j,0]<32.6):
            #do outermost edge
            iis[j] = mask_diff[0]
        
        elif len(mask_diff)==1:
            iis[j] = mask_diff[0]
        
        elif len(mask_diff)==0:
            iis[j] = iis[j-1]
        
        #save lon/lat for finding where to trim around TJRE
        lonshore[j] = lon_rho[j,iis[j]]
        latshore[j] = lat_rho[j,iis[j]]

    lon_diff = np.diff(lonshore)
    #cutoff before jumping at mouth of San Diego bay
    cutoff = np.argmax(np.abs(lon_diff))
    #find where mouth of TJRE is using the jump in longitude following shoreline
    TJRE_inds = np.where(np.abs(lon_diff[:cutoff])>0.0036)[0]
    TJ0 = TJRE_inds[0]
    TJ1 = TJRE_inds[-1]+1

    # rebuild iis and jjs around TJRE
    iis = np.concatenate([iis[:TJ0],iis[TJ1:cutoff]])
    jjs = np.concatenate([jjs[:TJ0],jjs[TJ1:cutoff]])

    # trim a little more - optional
    cutoff = 984
    iis = iis[:cutoff]
    jjs = jjs[:cutoff]

    return iis,jjs

def add_beaches_to_ax_RHS(ax,TJRE_line=False,xaxscale=0.02):
    beach_name_list = ['PB', 'TJRE', 'PTJ','SS','HdC','IB']
    beach_list = get_beach_location(beach_name_list)
    [x0,x1] = ax.get_xlim()
    [y0,y1] = ax.get_ylim()
    ax.set_ylim([2,y1])
    [y0,y1] = ax.get_ylim()
    # if type(x0)==np.float64:
    #     xaxscale = 0.1
    # else:
    #     xaxscale = 0.02
    for beach in beach_name_list:
        beach_r_km = 0.001*(beach_list[beach]['r']-beach_list['PB']['r'])
        if (beach_r_km>=y0) and (beach_r_km<=y1):
            if beach=='PB':
                beach_color='m'
                beach_marker = '<'
                text_color='k'
            elif beach=='TJRE':
                beach_color='yellowgreen'
                beach_marker = '<'
                text_color='k'
                if TJRE_line:
                    ax.axhline(y=0.001*(beach_list[beach]['r']-beach_list['PB']['r']),color=beach_color,linestyle='dashed',linewidth=1.0)
            else:
                beach_color='yellow'
                beach_marker = 'o'
                text_color='k'

            ax.plot(x1,0.001*(beach_list[beach]['r']-beach_list['PB']['r']),marker=beach_marker,mec='k',mfc=beach_color,markersize=8,clip_on=False,zorder=500)
            ax.set_xlim([x0,x1])
            
            ax.text(x1+xaxscale*(x1-x0),0.001*(beach_list[beach]['r']-beach_list['PB']['r']),beach,color=text_color,ha='left',va='center')

    

def wavebuoy_to_5m_isobath(Dbuoy):
    g = 9.81
    shorenormal = 263
    rho0 = 1028
    
    WaveBuoy = {}
    for var_name in ['Dwave','Lwave','Hwave']:
        #the first timestep has Hwave=0
        WaveBuoy[var_name] = Dbuoy[var_name][1:]
    
    WaveBuoy['Dprime'] = WaveBuoy['Dwave']-shorenormal
    WaveBuoy['Dprime_rad'] = WaveBuoy['Dprime']*np.pi/180
    WaveBuoy['k'] = (2*np.pi)/WaveBuoy['Lwave']
    WaveBuoy['h'] = Dbuoy['h0'][1:][Dbuoy['bj'],Dbuoy['bi']]+Dbuoy['zeta'][1:]
    kh = WaveBuoy['k']*WaveBuoy['h']
    WaveBuoy['cp'] = np.sqrt((g/WaveBuoy['k'])*np.tanh(kh))
    WaveBuoy['omega'] = WaveBuoy['cp']*WaveBuoy['k']
    WaveBuoy['Twave'] = 2*np.pi/WaveBuoy['omega']
    sechkh = 1/np.cosh(kh)
    WaveBuoy['cg'] = 0.5*(g*np.tanh(kh)+g*kh*sechkh**2)/np.sqrt(g*WaveBuoy['k']*np.tanh(kh))
    WaveBuoy['fit'] = 523.480

    Isobath5m = {}
    Isobath5m['h'] = 5
    Isobath5m['cp'] = np.sqrt(g*Isobath5m['h'])
    Isobath5m['Dprime_rad'] = np.arcsin(np.sin(WaveBuoy['Dprime']*np.pi/180)*(Isobath5m['cp']/WaveBuoy['cp']))
    Isobath5m['Dprime'] = Isobath5m['Dprime_rad']*180/np.pi
    Isobath5m['Dwave'] = Isobath5m['Dprime']+shorenormal
    Isobath5m['cg'] = np.sqrt(g*Isobath5m['h']) # in shallow water, cg -> root(gh), same as cp
    Isobath5m['Hwave'] = WaveBuoy['Hwave']*np.sqrt((WaveBuoy['cg']/Isobath5m['cg'])*(np.cos(WaveBuoy['Dprime_rad'])/np.cos(Isobath5m['Dprime_rad'])))
    Isobath5m['omega'] = WaveBuoy['omega']
    Isobath5m['k'] = Isobath5m['omega']/Isobath5m['cp']
    Isobath5m['Twave'] = 2*np.pi/Isobath5m['omega']
    Isobath5m['fit'] = 580.980
    
    for wave_loc in WaveBuoy,Isobath5m:
        wave_loc['E'] = -(1/16)*rho0 * g* wave_loc['Hwave']**2
        wave_loc['Sxy'] = wave_loc['E'] * (wave_loc['cg']/wave_loc['cp']) *  np.sin(wave_loc['Dprime_rad']) * np.cos(wave_loc['Dprime_rad'])
        denominator = wave_loc['fit']*wave_loc['Hwave']
        wave_loc['V'] = wave_loc['Sxy']/denominator
        wave_loc['linear_fit'] = 820.29
        
        # test Wright and Thompson 1983
        sigma = np.sqrt(g/wave_loc['h'])*wave_loc['Hwave']/4
        v_sigma = wave_loc['V']/sigma
        
    
    return WaveBuoy,Isobath5m

def gen_plot_props(fs_big=12,fs_small=10,lw_big=2,lw_small=1):
    plt.rc('xtick', labelsize=fs_small)
    plt.rc('ytick', labelsize=fs_small)
    plt.rc('xtick.major', size=8, pad=5, width=lw_small)
    plt.rc('ytick.major', size=8, pad=5, width=lw_small)
    plt.rc('axes', lw=lw_small)
    plt.rc('lines', lw=lw_big)
    plt.rc('font', size=fs_big)
    plt.rc('grid', color='g', ls='-', lw=lw_small, alpha=.3)
    plt.rc('axes', axisbelow=True)
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    
    #default is a 1/4 page figure. If 1/2 or full page is desired, x2
    # before using
    fw_mm = 95 # in millimeters
    fh_mm = 115 # in millimeters
    
    #conversion from millimeters to inches
    mm2in = 1/25.4
    fw = fw_mm*mm2in
    fh = fh_mm*mm2in
    
    return(fw,fh)

def find_nearest_ind_2D(x_mat,y_mat,x0,y0):
    # designed to help locate the i and j location of points on non-plaid 2D grids
    #
    # NOTE: pay attention to i,j vs j,i order
    # I used "i" to correspond to "x" and "j" to correspond to "y"
    # for ROMS grids in python, the correct indexing order will be [t,z,y,x]
    # So this code returns "i,j"
    # but the desired grid cell will be x = x_mat[j,i], y = y_mat[j,i]
    
    xdiff = x_mat-x0
    ydiff = y_mat-y0
    xydiff = np.sqrt(xdiff**2+ydiff**2)
    i = np.where(xydiff==xydiff.min())[1][0]
    j = np.where(xydiff==xydiff.min())[0][0]
    
    return(i,j)

def find_nearest_ind_2D(x_mat,y_mat,x0,y0):
    # designed to help locate the i and j location of points on non-plaid 2D grids
    #
    # NOTE: pay attention to i,j vs j,i order
    # I used "i" to correspond to "x" and "j" to correspond to "y"
    # for ROMS grids in python, the correct indexing order will be [t,z,y,x]
    # So this code returns "i,j"
    # but the desired grid cell will be x = x_mat[j,i], y = y_mat[j,i]
    
    xdiff = x_mat-x0
    ydiff = y_mat-y0
    xydiff = np.sqrt(xdiff**2+ydiff**2)
    i = np.where(xydiff==xydiff.min())[1][0]
    j = np.where(xydiff==xydiff.min())[0][0]
    
    return(i,j)

def rot_uv_to_shoreangle(u0,v0,shoreangle):
    #should support either scalar shoreangle or multidimensional shoreangle
    # as long as shoreangle is same shape as u, v
    # shoreangle should be in DEGREES, not radians
    
    #load in grid for rotating
    home = '/Users/elizabethbrasseale/Projects/Water quality/'
    ll_fn = home+'WQ_data/extractions2017/shoreline_variables_2017.p'
    D=pickle.load(open(ll_fn,'rb'))
    lat_rho = D['lat_rho'][:]
    lon_rho = D['lon_rho'][:]
    
    # bc shoreangle is calculated using lat/lon, first rotate velocity into N/S
    grid_angle = np.mean(np.arctan2(np.diff(lat_rho,axis=1),np.diff(lon_rho,axis=1)))
    w0 = u0 + 1j*v0
    w = w0 * np.exp(1j*grid_angle)

    #now rotate velocity along the shoreangles
    theta_SA = (shoreangle)*np.pi/180
    w = w * np.exp(1j*theta_SA)
    u = np.real(w)
    v = np.imag(w)
    
    return u,v

def willmott(m,o):
    """
    Calculates willmott skill score between two vectors, a model and a set of observations
    """
    # if len(m)!=len(o):
    #     error('Vectors must be the same length!');
    # end

    MSE = np.nanmean((m - o)**2)
    denom1 = abs(m - np.nanmean(o))
    denom2 = abs(o - np.nanmean(o))
    denom = np.nanmean((denom1 + denom2)**2)
    
    if denom==0:
        WS = 0
    else:
        WS = 1 - MSE/denom
    
    return WS

def smooth(vec,window_width):
    half_window = int(window_width/2)
    # cumsum_vec = np.cumsum(np.insert(vec, 0, 0))
    # vec_padded = np.pad(vec, (half_window, window_width-1-half_window), mode='edge')
    vec_padded = np.pad(vec, (half_window, window_width-1-half_window), mode='constant',constant_values=(np.mean(vec[:10]),np.mean(vec[-10:])))
    cumsum_vec = np.cumsum(np.insert(vec_padded,0,0))
    new_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return new_vec

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE

def ll2xy(lon, lat, lon0, lat0):
    """
    This converts lon, lat into meters relative to lon0, lat0.
    It should work for lon, lat scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y


def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

def filt_godin(data):
    """
    Input: 1D numpy array of HOURLY values
    Output: Array of the same size, filtered with 24-24-25 Godin filter,
        padded with nan's
    """
    filt = godin_shape()
    npad = np.floor(len(filt)/2).astype(int)
    smooth = np.convolve(data, filt, mode = 'same')
    smooth[:npad] = np.nan
    smooth[-npad:] = np.nan
    return smooth
    
def filt_godin_mat(data):
    """
    Input: ND numpy array of HOURLY, with time on axis 0.
    Output: Array of the same size, filtered with 24-24-25 Godin filter,
        padded with nan's
    """
    filt = godin_shape()
    n = np.floor(len(filt)/2).astype(int)
    sh = data.shape
    df = data.flatten('F')
    dfs = np.convolve(df, filt, mode = 'same')
    smooth = dfs.reshape(sh, order='F')
    smooth[:n,:] = np.nan
    smooth[-n:,:] = np.nan
    return smooth
    
def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def get_beach_location(beach_name_list):
    
    dir0 = '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/'
    fn = dir0 + 'extractions2017/shoreline_variables_2017.p'
    
    D = pickle.load(open(fn,'rb'))
    lonshore = D['lonshore'][:]
    latshore = D['latshore'][:]
    xshore = D['xshore'][:]
    yshore = D['yshore'][:]
    rshore = D['rshore'][:]
    lon_rho = D['lon_rho'][:]
    lat_rho = D['lat_rho'][:]
    
    # Identify some points of interest on the map
    PB = {'lat':32.446, 'name':'Punta Bandera'}
    TJRE = {'lat':32.553, 'name':'Tijuana River Estuary'}
    PTJ = {'lat':32.518, 'name':'Playas Tijuana'}
    IB = {'lat':32.579, 'name':'Imperial Beach Pier'}
    SS = {'lat':32.632, 'name':'Silver Strand State Beach'}
    HdC = {'lat':32.678, 'name':'Hotel del Coronado'}

    beach_list = {}

    for var in beach_name_list:
        beach_list[var] = locals()[var]

    for beach in beach_list:
        y_ind = np.argmin(np.abs(latshore-beach_list[beach]['lat']))
        beach_list[beach]['lon'] = lonshore[y_ind]
        beach_list[beach]['x'] = xshore[y_ind]
        beach_list[beach]['y'] = yshore[y_ind]
        beach_list[beach]['r'] = rshore[y_ind]
        
    return(beach_list)
   
def get_shoreline_models(model_name_list):

    dir0 = '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/'
    
    #load in data
    # #cside dye
    CSIDE_SAsm10 = {}
    cside_fn = dir0+'extractions2017/shoreline_dye_waves_05m_interp.p'
    Dcside = pickle.load(open(cside_fn,'rb'))
    CSIDE_SAsm10['dye'] = Dcside['dye_01'][:]
    CSIDE_SAsm10['y'] = Dcside['rshore'][:]

    #cside velocities
    cside_uv_fn = dir0 + 'extractions2017/shoreline_uv_05m_interp.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_SAsm10['v'] = Duv['v'][:]
    CSIDE_SAsm10['label']='CSIDE extracted data'# (rotated using\nsmoothed shoreangle (window=10))'

    CSIDE_PAsm100 = {}
    CSIDE_PAsm100['dye'] = Dcside['dye_01'][:,1:-1]
    CSIDE_PAsm100['y'] = Dcside['rshore'][:]
    cside_uv_fn = dir0 + 'extractions2017/shoreline_uv_05m_interp_sm100_PA.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_PAsm100['v'] = Duv['v'][:]
    CSIDE_PAsm100['label']='CSIDE extracted data (rotated using\nsmoothed principle axis (window=100))'

    CSIDE_SAsm100 = {}
    cside_uv_fn = dir0 + 'extractions2017/shoreline_dye_waves_uv_05m_interp_sm100_SA.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_SAsm100['dye'] = Duv['dye_01'][:,1:-1]
    CSIDE_SAsm100['y'] = Duv['rshore'][1:-1]
    CSIDE_SAsm100['v'] = Duv['v_rot_int'][:,1:-1]
    CSIDE_SAsm100['label']='3D COAWST model'
    
    CSIDE_2017_2018 = {}
    cside_uv_fn = dir0 + 'shoreline_variables_2017–2018.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_2017_2018['dye'] = Duv['dye_01_esg'][:,1:-1]
    CSIDE_2017_2018['y'] = Duv['rshore_esg'][1:-1]
    CSIDE_2017_2018['v'] = Duv['v_esg'][:,1:-1]
    CSIDE_2017_2018['label']='3D COAWST model'

    #adv-diff model wtih resolved alongshore-varying velocity w/ no smoothing
    AV_nosm = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_waves_postinterp.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_nosm['dye'] = da_model['c'][:]
    AV_nosm['y'] = da_model['y'][:]
    AV_nosm['v'] = da_model['v'][:]
    AV_nosm['label'] = 'Adv-Diff model with alongshore-varying input'#'\nusing smoothed shoreangle (window=10)'

    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_0km_tuned_subsampled.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_tuned['v'] = da_model['v'][:]
    AV_SAsm100_tuned['label'] = '1D Alongshore-varying wave advection model'#'\nusing smoothed shoreangle (window=10)'

    # #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    # AV_SAsm100_0km_tuned = {}
    # model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_0km_tuned_subsampled.p'
    # da_model = pickle.load(open(model_fn,'rb'))
    # AV_SAsm100_0km_tuned['dye'] = da_model['c'][:]
    # AV_SAsm100_0km_tuned['y'] = da_model['y'][1:-1]
    # AV_SAsm100_0km_tuned['v'] = da_model['v'][:]
    # AV_SAsm100_0km_tuned['label'] = 'Adv-Diff model with 0 km running average'#'\nusing smoothed shoreangle (window=100)'

    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_3km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_3km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_3km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_3km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_3km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_3km_tuned['label'] = 'Adv-Diff model with 3 km running average'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_4km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_4km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_4km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_4km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_4km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_4km_tuned['label'] = 'Adv-Diff model with 4 km running average'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_5km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_5km_tuned_subsampled.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_5km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_5km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_5km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_5km_tuned['label'] = '1D Alongshore-varying wave advection model (smoothed)'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_6km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_6km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_6km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_6km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_6km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_6km_tuned['label'] = 'Adv-Diff model with 6 km running average'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_7km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_7km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_7km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_7km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_7km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_7km_tuned['label'] = 'Adv-Diff model with 7 km running average'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_8km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_8km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_8km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_8km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_8km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_8km_tuned['label'] = 'Adv-Diff model with 8 km running average'#'\nusing smoothed shoreangle (window=100)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_10km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_10km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_10km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_10km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_10km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_10km_tuned['label'] = 'Adv-Diff model with 10 km running average'#'\nusing smoothed shoreangle (window=100)'

    #adv-diff model with CSIDE-extracted velocity input
    AV_recycled_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_recycled_input_tuned_subsampled.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_recycled_tuned['y'] = da_model['y'][1:-1]
    AV_recycled_tuned['dye']  = da_model['c'][:,1:-1]
    AV_recycled_tuned['v']  = da_model['v'][:,1:-1]
    AV_recycled_tuned['label'] = '1D COAWST advection model'

    #adv-diff model with 3km binned velocity
    AV_SAsm10_3kmbin = {}
    model_fn=dir0+'adv_diff_model/CSIDE_3kmbin_waves_postinterp.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_3kmbin['y'] = da_model['y'][:]
    AV_SAsm10_3kmbin['dye']  = da_model['c'][:]
    AV_SAsm10_3kmbin['v']  = da_model['v'][:]
    AV_SAsm10_3kmbin['label'] = 'Adv-Diff model with 3km-binned velocity input'

    #adv-diff model with 1km running average velocity
    AV_SAsm10_1kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_1km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_1kmavg['y'] = da_model['y'][:]
    AV_SAsm10_1kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_1kmavg['v']  = da_model['v'][:]
    AV_SAsm10_1kmavg['label'] = 'Adv-Diff model with 1km-running avg velocity input'

    #adv-diff model with 3km running average velocity
    AV_SAsm10_3kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_3km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_3kmavg['y'] = da_model['y'][:]
    AV_SAsm10_3kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_3kmavg['v']  = da_model['v'][:]
    AV_SAsm10_3kmavg['label'] = 'Adv-Diff model with 3km-running avg velocity input'

    #adv-diff model with 5km running average velocity
    AV_SAsm10_5kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_5km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_5kmavg['y'] = da_model['y'][:]
    AV_SAsm10_5kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_5kmavg['v']  = da_model['v'][:]
    AV_SAsm10_5kmavg['label'] = 'Adv-Diff model with 5km-running avg velocity input'


    #adv-diff model with uniform velocity
    U_tuned = {}
    model_fn2=dir0+'adv_diff_model/CSIDE_tuning_kd_-2.67E-05_notsubsampled.p'
    da_model2 = pickle.load(open(model_fn2,'rb'))
    U_tuned['y']= da_model2['y'][1:-1]
    U_tuned['dye'] = da_model2['c'][:,1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_tuned['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tuned['v'] = v0
    U_tuned['label'] = 'extracted near PB'
    
    #adv-diff model with uniform velocity
    U_tuned_2017_2018 = {}
    model_fn2=dir0+'adv_diff_model/CSIDE_tuning_kd_-2.67E-05_subsampled_2017–2018.p'
    da_model2 = pickle.load(open(model_fn2,'rb'))
    U_tuned_2017_2018['y']= da_model2['y'][1:-1]
    U_tuned_2017_2018['dye'] = da_model2['c'][:,1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_tuned_2017_2018['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tuned_2017_2018['v'] = v0
    U_tuned_2017_2018['label'] = '1D Wave advection model'
    
    #adv-diff model with uniform velocity
    U_tuned_ss = {}
    model_fn2=dir0+'adv_diff_model/CSIDE_tuning_kd_-2.67E-05_subsampled.p'
    da_model2 = pickle.load(open(model_fn2,'rb'))
    U_tuned_ss['y']= da_model2['y'][1:-1]
    U_tuned_ss['dye'] = da_model2['c'][:,1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_tuned_ss['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tuned_ss['v'] = v0
    U_tuned_ss['label'] = '1D Wave advection model'
    #Note: no performance difference from U_tuned

    #adv-diff model wtih resolved alongshore-varying velocity
    #but leaving waves uninterpolated, and interpolating v after calculating Sxy
    AV_PAsm100 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_PA.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_PAsm100['dye'] = da_model['c'][:]
    AV_PAsm100['y'] = da_model['y'][1:-1]
    AV_PAsm100['v'] = da_model['v'][:]
    AV_PAsm100['label'] = 'Adv-Diff model with alongshore-varying input\nusing smoothed PA (window=100)'

    #adv-diff model wtih resolved alongshore-varying velocity
    #but leaving waves uninterpolated, and interpolating v after calculating Sxy
    AV_SAsm100 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100['dye'] = da_model['c'][:]
    AV_SAsm100['y'] = da_model['y'][1:-1]
    AV_SAsm100['v'] = da_model['v'][:]
    AV_SAsm100['label'] = 'Adv-Diff model with alongshore-varying input\nusing smoothed shoreangles (window=100)'


    #alongshore diffusion model
    U_Kyy01 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy01_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy01['dye'] = da_model['c'][:,1:-1]
    U_Kyy01['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy01['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy01['v'] = v0
    U_Kyy01['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 1 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy03 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy03_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy03['dye'] = da_model['c'][:,1:-1]
    U_Kyy03['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy03['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy03['v'] = v0
    U_Kyy03['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 3 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy05 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy05_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy05['dye'] = da_model['c'][:,1:-1]
    U_Kyy05['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy05['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy05['v'] = v0
    U_Kyy05['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 5 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy07 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy07_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy07['dye'] = da_model['c'][:,1:-1]
    U_Kyy07['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy07['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy07['v'] = v0
    U_Kyy07['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 7 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy08 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy08_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy08['dye'] = da_model['c'][:,1:-1]
    U_Kyy08['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy08['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy08['v'] = v0
    U_Kyy08['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 8 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy09 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy09_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy09['dye'] = da_model['c'][:,1:-1]
    U_Kyy09['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy09['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy09['v'] = v0
    U_Kyy09['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 9 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy10 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy10_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy10['dye'] = da_model['c'][:,1:-1]
    U_Kyy10['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_Kyy10['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy10['v'] = v0
    U_Kyy10['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 10 m^{2}s^{-1}$'
    
    #shorenormal tuning model
    U_sa230 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal230.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa230['dye'] = da_model['c'][:,1:-1]
    U_sa230['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa230['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa230['v'] = v0
    U_sa230['label'] = r'Adv-Diff model shorenormal $230^{\circ}$'
    
    #shorenormal tuning model
    U_sa235 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal235.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa235['dye'] = da_model['c'][:,1:-1]
    U_sa235['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa235['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa235['v'] = v0
    U_sa235['label'] = r'Adv-Diff model shorenormal $235^{\circ}$'
    
    #shorenormal tuning model
    U_sa240 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal240.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa240['dye'] = da_model['c'][:,1:-1]
    U_sa240['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa240['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa240['v'] = v0
    U_sa240['label'] = r'Adv-Diff model shorenormal $240^{\circ}$'
    
    
    #shorenormal tuning model
    U_sa245 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal245.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa245['dye'] = da_model['c'][:,1:-1]
    U_sa245['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa245['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa245['v'] = v0
    U_sa245['label'] = r'Adv-Diff model shorenormal $245^{\circ}$'
    
    #shorenormal tuning model
    U_sa250 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal250.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa250['dye'] = da_model['c'][:,1:-1]
    U_sa250['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa250['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa250['v'] = v0
    U_sa250['label'] = r'Adv-Diff model shorenormal $250^{\circ}$'
    
    #shorenormal tuning model
    U_sa2525 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal252.5.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa2525['dye'] = da_model['c'][:,1:-1]
    U_sa2525['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa2525['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa2525['v'] = v0
    U_sa2525['label'] = r'Adv-Diff model shorenormal $252.5^{\circ}$'
    
    #shorenormal tuning model
    U_sa255 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal255.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa255['dye'] = da_model['c'][:,1:-1]
    U_sa255['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa255['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa255['v'] = v0
    U_sa255['label'] = r'Adv-Diff model shorenormal $255^{\circ}$'
    
    #shorenormal tuning model
    U_sa2575 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal257.5.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa2575['dye'] = da_model['c'][:,1:-1]
    U_sa2575['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa2575['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa2575['v'] = v0
    U_sa2575['label'] = r'Adv-Diff model shorenormal $257.5^{\circ}$'
    
    #shorenormal tuning model
    U_sa260 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal260.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa260['dye'] = da_model['c'][:,1:-1]
    U_sa260['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa260['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa260['v'] = v0
    U_sa260['label'] = r'Adv-Diff model shorenormal $260^{\circ}$'

    #shorenormal tuning model
    U_sa2625 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal262.5.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa2625['dye'] = da_model['c'][:,1:-1]
    U_sa2625['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa2625['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa2625['v'] = v0
    # U_sa2625['label'] = r'Adv-Diff model shorenormal $262.5^{\circ}$'
    U_sa2625['label'] = r'Extracted near PB'
    
    #shorenormal tuning model
    U_sa265 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal265.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa265['dye'] = da_model['c'][:,1:-1]
    U_sa265['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa265['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa265['v'] = v0
    U_sa265['label'] = r'Adv-Diff model shorenormal $265^{\circ}$'

    #shorenormal tuning model
    U_sa2675 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal267.5.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa2675['dye'] = da_model['c'][:,1:-1]
    U_sa2675['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa2675['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa2675['v'] = v0
    U_sa2675['label'] = r'Adv-Diff model shorenormal $267.5^{\circ}$'
    
    #shorenormal tuning model
    U_sa270 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal270.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_sa270['dye'] = da_model['c'][:,1:-1]
    U_sa270['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_sa270['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_sa270['v'] = v0
    U_sa270['label'] = r'Adv-Diff model shorenormal $270^{\circ}$'
    
    # tuned bottom drag
    U_tunedBD = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_Kyy7.0_kd-2.67E-05_shorenormal255_bottomdragtuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedBD['dye'] = da_model['c'][:,1:-1]
    U_tunedBD['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedBD['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedBD['v'] = v0
    U_tunedBD['label'] = r'Uniform V, tuned bottom drag coefficient'
    
    # tuned linear drag coefficient & shorenormal
    U_tunedLDC_SN = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_IBP_kd-2.67E-05_tuned_shorenormal_lineardragcoef.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedLDC_SN['dye'] = da_model['c'][:,1:-1]
    U_tunedLDC_SN['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedLDC_SN['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedLDC_SN['v'] = v0
    U_tunedLDC_SN['label'] = 'IBP bouy'
    
    # tuned linear drag coefficient & shorenormal
    U_tunedLDC_SN_Lsz = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_IBP_kd-2.67E-05_tuned_shorenormal_lineardragcoef_varLsz.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedLDC_SN_Lsz['dye'] = da_model['c'][:,1:-1]
    U_tunedLDC_SN_Lsz['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedLDC_SN_Lsz['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedLDC_SN_Lsz['v'] = v0
    U_tunedLDC_SN_Lsz['label'] = 'IBP bouy'
    
    # tuned linear drag coefficient & shorenormal
    U_tunedLDC_SN_Lsz_kd = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_IBP_kd-1.40E-05_tuned_shorenormal_lineardragcoef_varLsz.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedLDC_SN_Lsz_kd['dye'] = da_model['c'][:,1:-1]
    U_tunedLDC_SN_Lsz_kd['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedLDC_SN_Lsz_kd['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedLDC_SN_Lsz_kd['v'] = v0
    U_tunedLDC_SN_Lsz_kd['label'] = 'IBP bouy, tuned kd'
    
    # tuned linear drag coefficient & shorenormal
    U_tunedLDC_SN_Lsz_kd_lowPB = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_IBP_kd-1.40E-05_tuned_shorenormal_lineardragcoef_varLsz_lowPB_in.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedLDC_SN_Lsz_kd_lowPB['dye'] = da_model['c'][:,1:-1]
    U_tunedLDC_SN_Lsz_kd_lowPB['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedLDC_SN_Lsz_kd_lowPB['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedLDC_SN_Lsz_kd_lowPB['v'] = v0
    U_tunedLDC_SN_Lsz_kd_lowPB['label'] = 'IBP bouy, tuned kd, PB*1/10'
    
    # tuned linear drag coefficient & shorenormal
    U_tunedLDC_SN_Lsz_kd_medPB = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_IBP_kd-1.40E-05_tuned_shorenormal_lineardragcoef_varLsz_medPB_in.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_tunedLDC_SN_Lsz_kd_medPB['dye'] = da_model['c'][:,1:-1]
    U_tunedLDC_SN_Lsz_kd_medPB['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_tunedLDC_SN_Lsz_kd_medPB['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tunedLDC_SN_Lsz_kd_medPB['v'] = v0
    U_tunedLDC_SN_Lsz_kd_medPB['label'] = 'IBP bouy, tuned kd, PB*1/5'
    
    # tuned linear drag coefficient & shorenormal
    U_isobath5m = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_5miso_kd-1.40E-05_tuned_shorenormal_medPB_in.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m['dye'] = da_model['c'][:,1:-1]
    U_isobath5m['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_isobath5m['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m['v'] = v0
    U_isobath5m['label'] = 'IBP bouy, propagated to 5m isobath'
    
    # tuned  shorenormal w/ tuned C0
    U_isobath5m_C0 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_5miso_kd-1.40E-05_tuned_shorenormal_PB_in0.012.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m_C0['dye'] = da_model['c'][:,1:-1]
    U_isobath5m_C0['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_isobath5m_C0['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m_C0['v'] = v0
    U_isobath5m_C0['label'] = 'IBP bouy, propagated to 5m isobath'
    
    # tuned linear drag coefficient & shorenormal
    U_isobath5m_r = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_5miso_rayleigh_kd-1.40E-05_tuned_shorenormal_medPB_in.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m_r['dye'] = da_model['c'][:,1:-1]
    U_isobath5m_r['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_isobath5m_r['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m_r['v'] = v0
    U_isobath5m_r['label'] = 'IBP bouy, propagated to 5m isobath, Rayleigh friction'
    
    # tuned linear drag coefficient & shorenormal w/ tuned C0
    U_isobath5m_r_C0 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_uniformV_5miso_rayleigh_kd-1.40E-05_tuned_shorenormal_PB_in0.02.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m_r_C0['dye'] = da_model['c'][:,1:-1]
    U_isobath5m_r_C0['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_isobath5m_r_C0['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m_r_C0['v'] = v0
    U_isobath5m_r_C0['label'] = 'IBP bouy, propagated to 5m isobath, Rayleigh friction'
    
    # autotuned sawc drag coefficient & shorenormal w/ tuned C0
    U_isobath5m_sawc_autotune_kd = {}
    model_fn=dir0+'adv_diff_model/autotuned_sawc_kd1.3000E-05_uniformv_5miso_PB_in0.008.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m_sawc_autotune_kd['dye'] = da_model['c'][:,1:-1]
    U_isobath5m_sawc_autotune_kd['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj = U_isobath5m_sawc_autotune_kd['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m_sawc_autotune_kd['v'] = v0
    U_isobath5m_sawc_autotune_kd['label'] = 'IBP bouy, propagated to 5m isobath, SAWV, autotuned kd=-5.5e-6'
    
    # autotuned rayleigh drag coefficient & shorenormal w/ tuned C0
    U_isobath5m_R_autotune_kd = {}
    model_fn=dir0+'adv_diff_model/autotuned_rayleigh_kd1.3000E-05_uniformv_5miso_PB_in0.010.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_isobath5m_R_autotune_kd['dye'] = da_model['c'][:,1:-1]
    U_isobath5m_R_autotune_kd['y'] = da_model['y'][1:-1]
    v00 = da_model['v'][:]
    nt,nj =U_isobath5m_R_autotune_kd['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_isobath5m_R_autotune_kd['v'] = v0
    U_isobath5m_R_autotune_kd['label'] = 'IBP bouy, propagated to 5m isobath, Rayleigh friction, autotuned kd=-6.5e-6'
    
    # tuned linear drag coefficient & shorenormal
    AV_recycled_tuned_medPB = {}
    model_fn=dir0+'adv_diff_model/CSIDE_recycled_input_tuned_subsampled_medPB_in.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_recycled_tuned_medPB['dye'] = da_model['c'][:,1:-1]
    AV_recycled_tuned_medPB['y'] = da_model['y'][1:-1]
    AV_recycled_tuned_medPB['v'] = da_model['v'][:,1:-1]
    AV_recycled_tuned_medPB['label'] = 'alongshore-varying recycled SD Bight v'
    
    # tuned linear drag coefficient & shorenormal
    AV_recycled_tuned_C0 = {}
    model_fn=dir0+'adv_diff_model/autotuned_recycled_kd1.3000E-05_AVv_5miso_PB_in0.011.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_recycled_tuned_C0['dye'] = da_model['c'][:,1:-1]
    AV_recycled_tuned_C0['y'] = da_model['y'][1:-1]
    AV_recycled_tuned_C0['v'] = da_model['v'][:,1:-1]
    AV_recycled_tuned_C0['label'] = 'alongshore-varying recycled SD Bight v'

    model_dict_list = {}

    for var in model_name_list:
        model_dict_list[var] = locals()[var]
    
    return(model_dict_list)