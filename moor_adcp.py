import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
from scipy import interpolate
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import seawater as sw

initials = ['NB'] #adcp is at North Buoy/Hansville

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_0 = '/Users/erinbroatch/Documents/GitHub/LO_output/extract/cas6_v3_lo8b/moor/orca_adcp_eb/'
fn = {'NB':fn_0 + 'NB_2021.01.01_2021.12.31.nc'}

fn_out_0 = '/Users/erinbroatch/Documents/Research/Data/LO_orca_moor/datasets/'
fn_out = {}
#can add more versions here
for j in initials:
    fn_out[j] = fn_out_0 + j + '_adcp_ds.nc'

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
    #press=np.linspace(plim[j][0], plim[j][1], plim[j][1]-plim[j][0]+1)
    press=np.linspace(plim[j][0], plim[j][1], 40)

    #make secondary depth axis
    depth_sw = sw.dpth(press, lonlat[j][1])

    #load data
    raw_ds=xr.open_dataset(fn[j])

    #make time index
    Btime = raw_ds.ocean_time.values
    dti = pd.to_datetime(Btime) #make datetimeindex

    #get depth array referenced to surface
    Bdepth = raw_ds.z_rho.values
    Bdepth_w = raw_ds.z_w.values
    Bsurf = raw_ds.z_w.values[:,[30]]
    depth_surf = Bdepth-Bsurf
    depth_surf_w = Bdepth_w - Bsurf

    #get variables
    u_raw = raw_ds.u.values
    v_raw = raw_ds.v.values
    w_raw = raw_ds.w.values

    #initialize arrays for interpolated data
    u_interp = np.zeros((len(dti),len(depth_sw)))
    v_interp = np.zeros((len(dti),len(depth_sw)))
    w_interp = np.zeros((len(dti),len(depth_sw)))

    #interpolate each profile in 1D onto regular grid matching orca data
    #need to flip input arrays so that depths are positive and increasing
    for k in range(len(dti)):
        u_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(u_raw[k,:]),right=np.nan)
        v_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf[k,:]),np.flip(v_raw[k,:]),right=np.nan)
        w_interp[k,:] = np.interp(depth_sw,np.flip(-depth_surf_w[k,:]),np.flip(w_raw[k,:]),right=np.nan)

    #Make DataArrays
    u = xr.DataArray(data=u_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='u velocity',units='m/s'))

    v = xr.DataArray(data=v_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='v velocity',units='m/s'))

    w = xr.DataArray(data=w_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='w velocity',units='m/s'))

    #Make Dataset
    ds = xr.Dataset()
    ds['u']=u
    ds['v']=v
    ds['w']=w
    ds=ds.assign_coords({'pressure':('z', press)}) #add pressure as alternate z-coordinate
    ds.pressure.attrs['long_name']='Pressure' #add attributes to dataset coordinates
    ds.pressure.attrs['units']='dbar'
    ds.time.attrs['long_name']='Time'
    ds.depth.attrs['long_name']='Depth'
    ds.depth.attrs['units']='m'

    # # Export as netCDF files !!!!!should change this to have a with open statement!!!!!
    ds.to_netcdf(fn_out[j])