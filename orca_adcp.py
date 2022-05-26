import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
from scipy import interpolate
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import seawater as sw

initials = ['NB']

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_0 = '/Users/erinbroatch/Documents/Research/Data/orca_adcp/'
fn = {'NB':fn_0 + 'HORCA_ADCP_2021.nc'}

fn_out_0 = '/Users/erinbroatch/Documents/Research/Data/orca_adcp/datasets/'
fn_out = {}
fn_out_hourly = {}
for j in initials:
    fn_out[j] = fn_out_0 + j + '_adcp_ds.nc'
    fn_out_hourly[j] = fn_out_0 + j + '_adcp_ds_hourly.nc'

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

    #time coordinnate
    dtnum_raw = raw_ds.dtnum
    dtnum=dtnum_raw.where(dtnum_raw<1e6).values #get rid of big dtnums

    #depth coordinate
    bin_depth = raw_ds.bin_depth.values

    #velocities
    u_raw = raw_ds.u_vel.values
    v_raw = raw_ds.v_vel.values
    w_raw = raw_ds.w_vel.values
    echo_raw = raw_ds.ecAve.values

    #delete data where time coordinate is NaN
    badtimes=np.where(np.isnan(dtnum))
    dtnum = np.delete(dtnum, badtimes)
    u_raw = np.delete(u_raw, badtimes, axis=0)
    v_raw = np.delete(v_raw, badtimes, axis=0)
    w_raw = np.delete(w_raw, badtimes, axis=0)
    echo_raw = np.delete(echo_raw, badtimes, axis=0)

    #make datetimeindex
    dtnum2000 = 730486 #reference point (Jan 1 2000) for converting to datetime
    dtnumref2000=dtnum-dtnum2000
    Bdatetime = []
    for i in range(len(dtnum)):
        Bdatetime.append(dt.timedelta(days=dtnumref2000[i])+dt.datetime(2000,1,1))
    dti = pd.to_datetime(Bdatetime).round('1s')

    #initialize arrays for interpolated data
    u_interp = np.zeros((len(dti),len(depth_sw)))
    v_interp = np.zeros((len(dti),len(depth_sw)))
    w_interp = np.zeros((len(dti),len(depth_sw)))
    echo_interp = np.zeros((len(dti),len(depth_sw)))

    #interpolate each profile in 1D onto regular grid matching orca data
    for k in range(len(dti)):
        u_interp[k,:] = np.interp(depth_sw,bin_depth,u_raw[k,:],right=np.nan)
        v_interp[k,:] = np.interp(depth_sw,bin_depth,v_raw[k,:],right=np.nan)
        w_interp[k,:] = np.interp(depth_sw,bin_depth,w_raw[k,:],right=np.nan)
        echo_interp[k,:] = np.interp(depth_sw,bin_depth,echo_raw[k,:],right=np.nan)

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

    echo = xr.DataArray(data=echo_interp,
        coords={'time':('time',dti), 'depth':('z', depth_sw)},
        dims=['time', 'z'],
        attrs=dict(long_name='Average echo intensity',units='counts'))

    #Make Dataset
    ds = xr.Dataset()
    ds['u']=u
    ds['v']=v
    ds['w']=w
    ds['echo']=echo
    ds=ds.assign_coords({'pressure':('z', press)}) #add pressure as alternate z-coordinate
    ds.pressure.attrs['long_name']='Pressure' #add attributes to dataset coordinates
    ds.pressure.attrs['units']='dbar'
    ds.time.attrs['long_name']='Time'
    ds.depth.attrs['long_name']='Depth'
    ds.depth.attrs['units']='m'

    if j=='NB':
        tdrop=ds.time.sel(time=slice("2021-05-24", "2021-07-17"))
        ds=ds.drop_sel(time=tdrop)

    # Make an hourly dataset
    start_date=str(ds.time[0].values)
    end_date=str(ds.time[-1].values)
    hourly_index=pd.date_range(start_date,end_date,freq='H')
    ds_hourly=ds.reindex({'time': hourly_index})
    ds_hourly=ds_hourly.transpose('z', 'time')

    # # Export as netCDF files !!!!!should change this to have a with open statement!!!!!
    ds.to_netcdf(fn_out[j])
    ds_hourly.to_netcdf(fn_out_hourly[j])