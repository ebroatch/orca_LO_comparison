import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
from scipy import interpolate
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import seawater as sw

initials = ['CI', 'PW', 'NB', 'DB', 'HP', 'TW'] #ADD NB DATA

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_0 = '/Users/erinbroatch/Documents/GitHub/LO_output/extract/cas6_v3_lo8b/moor/orca_eb/'
fn = {'CI':fn_0 + 'CI_2017.01.01_2021.12.31.nc',
    'PW':fn_0 + 'PW_2017.01.01_2021.12.31.nc',
    'NB':fn_0 + 'NB_2017.01.01_2021.12.31.nc',
    'DB':fn_0 + 'DB_2017.01.01_2021.12.31.nc',
    'HP':fn_0 + 'HP_2017.01.01_2021.12.31.nc',
    'TW':fn_0 + 'TW_2017.01.01_2021.12.31.nc'}

fn_out_0 = '/Users/erinbroatch/Documents/Research/Data/LO_orca_moor/datasets/'
fn_out = {}
#can add more versions here
for j in initials:
    fn_out[j] = fn_out_0 + j + '_ds.nc'

plim = {'CI':[3, 46],
    'PW':[3, 91],
    'NB':[3, 101],
    'DB':[3, 106],
    'HP':[3, 121],
    'TW':[3, 31]}

lonlat = {'CI': [-122.7300, 47.2800],
    'PW': [-122.3972, 47.7612],
    'NB': [-122.6270, 47.9073],
    'DB': [-122.8029, 47.8034],
    'HP': [-123.1126, 47.4218],
    'TW': [-123.0083, 47.3750]}

for j in initials:
    #make pressure axis
    press=np.linspace(plim[j][0], plim[j][1], plim[j][1]-plim[j][0]+1)

    #make secondary depth axis
    depth_sw = sw.dpth(press, lonlat[j][1])

    #load data
    raw_ds=xr.open_dataset(fn[j])

    #make time index
    Btime = raw_ds.ocean_time.values
    dti = pd.to_datetime(pd.to_datetime(Btime).date) #strip times and make datetimeindex

    #get depth array referenced to surface
    Bdepth = raw_ds.z_rho.values
    Bsurf = raw_ds.z_w.values[:,[30]]
    depth_surf = Bdepth-Bsurf

    #get variables
    Bsal = raw_ds.salt.values
    Btemp = raw_ds.temp.values
    Boxy = raw_ds.oxygen.values #oxygen in uM
    Boxy_mgL = Boxy * 32/1000 #convert to mg/L
    Bnitrate = raw_ds.NO3.values #nitrate in uM

    #initialize arrays for interpolated data
    Bsal_interp = np.zeros((len(dti),len(depth_sw)))
    Btemp_interp = np.zeros((len(dti),len(depth_sw)))
    Boxy_mgL_interp = np.zeros((len(dti),len(depth_sw)))
    Bnitrate_interp = np.zeros((len(dti),len(depth_sw)))

    #interpolate each profile in 1D onto regular grid matching orca data
    #need to flip input arrays so that depths are positive and increasing
    for k in range(len(dti)):
        Bsal_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(Bsal[k,:]),right=np.nan)
        Btemp_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(Btemp[k,:]),right=np.nan)
        Boxy_mgL_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(Boxy_mgL[k,:]),right=np.nan)
        Bnitrate_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(Bnitrate[k,:]),right=np.nan)

    #Make DataArrays
    sal = xr.DataArray(data=Bsal_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='Salinity',units='psu'))

    temp = xr.DataArray(data=Btemp_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='Temperature',units='degC'))

    oxy = xr.DataArray(data=Boxy_mgL_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='Oxygen',units='mg/L'))

    nitrate = xr.DataArray(data=Bnitrate_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='Nitrate',units='umol'))

    #Make Dataset
    ds = xr.Dataset()
    ds['sal']=sal
    ds['temp']=temp
    ds['oxy']=oxy
    ds['nitrate']=nitrate
    ds=ds.assign_coords({'pressure':('z', press)}) #add pressure as alternate z-coordinate
    ds.pressure.attrs['long_name']='Pressure' #add attributes to dataset coordinates
    ds.pressure.attrs['units']='dbar'
    ds.time.attrs['long_name']='Time'
    ds.depth.attrs['long_name']='Depth'
    ds.depth.attrs['units']='m'

    # # Export as netCDF files !!!!!should change this to have a with open statement!!!!!
    ds.to_netcdf(fn_out[j])