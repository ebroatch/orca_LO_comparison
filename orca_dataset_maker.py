import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import seawater as sw

initials = ['CI', 'PW', 'NB', 'DB', 'HP', 'TW']

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_0 = '/Users/erinbroatch/Documents/Research/Data/ORCA/'
fn = {'CI':fn_0 + 'CI_CTD_data_bin_web.mat',
    'PW':fn_0 + 'PW_CTD_data_bin_web.mat',
    'NB':fn_0 + 'HC_NB_CTD_data_bin_web.mat',
    'DB':fn_0 + 'HC_DB_CTD_data_bin_web.mat',
    'HP':fn_0 + 'HC_HP_CTD_data_bin_web.mat',
    'TW':fn_0 + 'HC_TW_CTD_data_bin_web.mat'}

fn_out_0 = '/Users/erinbroatch/Documents/Research/Data/ORCA/datasets/'
fn_out = {}
fn_out_daily = {}
fn_out_clim = {}
fn_out_seas = {}
for j in initials:
    fn_out[j] = fn_out_0 + j + '_ds.nc'
    fn_out_daily[j] = fn_out_0 + j + '_ds_daily.nc'
    fn_out_clim[j] = fn_out_0 + j + '_ds_clim.nc'
    fn_out_seas[j] = fn_out_0 + j + '_ds_seas.nc'

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
    mat = scipy.io.loadmat(fn[j])
    # if j=='CI':
    Btimefull = mat['Btime'][1:,:] #cut off top row
    indtime = np.argmax(np.isfinite(Btimefull),axis=0) #find index of first non-nan in column
    Btime = np.asarray([Btimefull[item, enum] for enum,item in enumerate(indtime)]) #use first non-nan (if available) as time coord
    #Btime = mat['Btime'][1,:] #second row of time array to use as time coordinate
    Bsal=mat['Bsal'][1:,:] #cut off top row
    Btemp=mat['Btemp'][1:,:]
    Boxy_mgL_cal=mat['Boxy_mgL_cal'][1:,:]
    Bnitrate_cal=mat['Bnitrate_cal'][1:,:]
    Bfluor_cal=mat['Bfluor_cal'][1:,:]
    Bpar=mat['Bpar'][1:,:]
    Bdepth=mat['Bdepth'][1:,:]
    # else:
    #     Btime = mat['Btime'][0,:]
    #     Bsal=mat['Bsal']
    #     Btemp=mat['Btemp']
    #     Boxy_mgL_cal=mat['Boxy_mgL_cal']
    #     Bnitrate_cal=mat['Bnitrate_cal']
    #     Bfluor_cal=mat['Bfluor_cal']
    #     Bpar=mat['Bpar']

    #delete columns where time axis is NaN
    badtimes=np.where(np.isnan(Btime))
    Btime = np.delete(Btime, badtimes)
    Bsal = np.delete(Bsal, badtimes, axis=1)
    Btemp = np.delete(Btemp, badtimes, axis=1)
    Boxy_mgL_cal = np.delete(Boxy_mgL_cal, badtimes, axis=1)
    Bnitrate_cal = np.delete(Bnitrate_cal, badtimes, axis=1)
    Bfluor_cal = np.delete(Bfluor_cal, badtimes, axis=1)
    Bpar = np.delete(Bpar, badtimes, axis=1)
    Bdepth = np.delete(Bdepth, badtimes, axis=1)

    # Convert decimal days to datetimeindex
    offset=dt.timedelta(hours=-8)
    tz_PDT=dt.timezone(offset, name="PDT")
    Bdatetime_utc=[]
    for i in range(len(Btime)):
        Bdatetime=dt.timedelta(days=Btime[i])+dt.datetime(2000,1,1,tzinfo=tz_PDT)
        Bdatetime_utc.append(Bdatetime.astimezone(dt.timezone.utc).replace(tzinfo=None))
    dti = pd.to_datetime(Bdatetime_utc)

    #Make DataArrays

    sal = xr.DataArray(data=Bsal,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Salinity',units='psu'))
    sal=sal.where(sal<36) #reasonable limits

    temp = xr.DataArray(data=Btemp,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Temperature',units='degC'))
    temp=temp.where(temp>0) #reasonable limits

    oxy = xr.DataArray(data=Boxy_mgL_cal,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Oxygen',units='mg/L',description='Calibrated oxygen'))
    oxy=oxy.where(oxy>0) #can't be negative

    nitrate = xr.DataArray(data=Bnitrate_cal,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Nitrate',units='umol',description='Calibrated nitrate'))
    nitrate=nitrate.where(nitrate>0) #can't be negative

    fluor = xr.DataArray(data=Bfluor_cal,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Fluorometer',units='mg/m^3',description='Calibrated fluorometer'))
    fluor=fluor.where(fluor!=-555) #error code
    fluor=fluor.where(fluor>0) #can't be negative

    par = xr.DataArray(data=Bpar,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Photosynthetically active radiation',units='uEinstein/m^2s'))
    par=par.where(par!=-555) #error code
    par=par.where(par>0) #can't be negative

    depth_obs = xr.DataArray(data=Bdepth,
        coords={'pressure':('z', press), 'time':('time',dti)},
        dims=['z','time'],
        attrs=dict(long_name='Observed depth',units='m'))

    #Make Dataset
    ds = xr.Dataset()
    ds['sal']=sal
    ds['temp']=temp
    ds['oxy']=oxy
    ds['nitrate']=nitrate
    ds['fluor']=fluor
    ds['par']=par
    ds['depth_obs']=depth_obs
    ds=ds.assign_coords({'depth':('z', depth_sw)}) #add depth as alternate z-coordinate
    ds.pressure.attrs['long_name']='Pressure' #add attributes to dataset coordinates
    ds.pressure.attrs['units']='dbar'
    ds.time.attrs['long_name']='Time'
    ds.depth.attrs['long_name']='Depth'
    ds.depth.attrs['units']='m'

    # Make a daily dataset
    ds_daily=ds.groupby('time.date').mean('time',keep_attrs=True)
    start_date=str(ds_daily.date[0].values) #make time coordinate a datetime again
    end_date=str(ds_daily.date[-1].values)
    ds_daily=ds_daily.assign_coords(date=pd.to_datetime(ds_daily.date))
    daily_index=pd.date_range(start_date,end_date,freq='D')
    ds_daily=ds_daily.reindex({'date': daily_index})
    ds_daily=ds_daily.rename({'date': 'time'})
    ds.time.attrs['long_name']='Time' #change this??
    ds_daily=ds_daily.transpose('z', 'time')

    # Make a climatology
    ds_clim=ds_daily.groupby('time.dayofyear').mean('time',keep_attrs=True)
    ds_clim=ds_clim.rename({'dayofyear': 'yearday'})
    ds_clim.yearday.attrs['long_name']='Yearday'
    ds_clim=ds_clim.transpose('z', 'yearday')

    # Make seasonal averages #should maybe change weightings here
    ds_seas=ds_daily.groupby('time.season').mean('time',keep_attrs=True)
    ds_seas.season.attrs['long_name']='Season'
    ds_seas=ds_seas.transpose('z', 'season')

    # # Export as netCDF files !!!!!should change this to have a with open statement!!!!!
    ds.to_netcdf(fn_out[j])
    ds_daily.to_netcdf(fn_out_daily[j])
    ds_clim.to_netcdf(fn_out_clim[j])
    ds_seas.to_netcdf(fn_out_seas[j])