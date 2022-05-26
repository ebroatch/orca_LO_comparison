import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
#import cmocean.cm as cmo
import datetime as dt
import netCDF4 as nc

plt.close('all')

initials = ['CI', 'PW', 'NB', 'DB', 'HP', 'TW']

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

# ORCA data
fn_0 = '/Users/erinbroatch/Documents/Research/Data/ORCA/datasets/'
fn_daily = {}
for j in initials:
    fn_daily[j] = fn_0 + j + '_ds_daily.nc'
ds_daily = {}
for j in initials:
    ds_daily[j] = xr.open_dataset(fn_daily[j])

# Interpolated LO data
fn_lo_0 = '/Users/erinbroatch/Documents/Research/Data/LO_orca_moor/datasets/'
fn_lo = {}
#can add more versions here
for j in initials:
    fn_lo[j] = fn_lo_0 + j + '_ds.nc'
ds_lo = {}
for j in initials:
    ds_lo[j] = xr.open_dataset(fn_lo[j])

# Location to save figures
fn_out = '/Users/erinbroatch/Documents/Research/Figures/2022_January_10/orca_lo_lineplots/'

for j in initials:

    #Find the "bottom" and size of 20% depth
    depthindseries=ds_daily[j]['sal'].notnull().sum(dim='z')
    depthindmean=depthindseries.where(depthindseries!=0).mean()
    ind20=int(np.round(depthindmean.values/5))
    depthind=int(np.round(depthindmean))

    #Slice top and bottom 20%
    saltop_orca=ds_daily[j]['sal'].isel(z=slice(0,ind20)).mean(dim='z')
    saltop_orca.attrs['long_name']='Salinity'
    saltop_orca.attrs['units']='psu'
    saltop_orca.attrs['Description']='ORCA mean salinity in top 20 percent of water column'
    saltop_lo=ds_lo[j]['sal'].isel(z=slice(0,ind20)).mean(dim='z')
    saltop_lo.attrs['long_name']='Salinity'
    saltop_lo.attrs['units']='psu'
    saltop_lo.attrs['Description']='LO mean salinity in top 20 percent of water column'
    salbot_orca=ds_daily[j]['sal'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')
    salbot_orca.attrs['long_name']='Salinity'
    salbot_orca.attrs['units']='psu'
    salbot_orca.attrs['Description']='ORCA mean salinity in bottom 20 percent of water column'
    salbot_lo=ds_lo[j]['sal'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')
    saltop_orca.attrs['long_name']='Salinity'
    salbot_lo.attrs['units']='psu'
    salbot_lo.attrs['Description']='LO mean salinity in bottom 20 percent of water column'

    temptop_orca=ds_daily[j]['temp'].isel(z=slice(0,ind20)).mean(dim='z')
    temptop_orca.attrs['long_name']='Temperature'
    temptop_orca.attrs['units']='degC'
    temptop_orca.attrs['Description']='ORCA mean temperature in top 20 percent of water column'
    temptop_lo=ds_lo[j]['temp'].isel(z=slice(0,ind20)).mean(dim='z')
    temptop_lo.attrs['long_name']='Temperature'
    temptop_lo.attrs['units']='degC'
    temptop_lo.attrs['Description']='LO mean temperature in top 20 percent of water column'
    tempbot_orca=ds_daily[j]['temp'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')
    tempbot_orca.attrs['long_name']='Temperature'
    tempbot_orca.attrs['units']='degC'
    tempbot_orca.attrs['Description']='ORCA mean temperature in bottom 20 percent of water column'
    tempbot_lo=ds_lo[j]['temp'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')
    tempbot_lo.attrs['long_name']='Temperature'
    tempbot_lo.attrs['units']='degC'
    tempbot_lo.attrs['Description']='LO mean temperature in bottom 20 percent of water column'

    oxytop_orca=ds_daily[j]['oxy'].isel(z=slice(0,ind20)).mean(dim='z')
    oxytop_lo=ds_lo[j]['oxy'].isel(z=slice(0,ind20)).mean(dim='z')
    oxybot_orca=ds_daily[j]['oxy'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')
    oxybot_lo=ds_lo[j]['oxy'].isel(z=slice(depthind-ind20,depthind)).mean(dim='z')

    #Make plot
    fig, axs = plt.subplots(3,1,sharex=True)

    #Salinity
    saltop_lo.plot(lw=2,color='paleturquoise',label='LO top 20%',ax=axs[0])
    salbot_lo.plot(lw=2,color='silver',label='LO bottom 20%',ax=axs[0])
    saltop_orca.plot(lw=1,color='teal',label='ORCA top 20%',ax=axs[0])
    salbot_orca.plot(lw=1,color='k',label='ORCA bottom 20%',ax=axs[0])
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_ylabel('Salinity [psu]')
    axs[0].set_xlabel('')
    axs[0].set_xlim([dt.datetime(2017,1,1),dt.datetime(2021,12,31)])

    #Temperature
    temptop_lo.plot(lw=2,color='paleturquoise',label='LO top 20%',ax=axs[1])
    tempbot_lo.plot(lw=2,color='silver',label='LO bottom 20%',ax=axs[1])
    temptop_orca.plot(lw=1,color='teal',label='ORCA top 20%',ax=axs[1])
    tempbot_orca.plot(lw=1,color='k',label='ORCA bottom 20%',ax=axs[1])
    #axs[1].legend()
    axs[1].grid(True)
    axs[1].set_ylabel('Temperature [degC]')
    axs[1].set_xlabel('')

    #Oxygen
    oxytop_lo.plot(lw=2,color='paleturquoise',label='LO top 20%',ax=axs[2])
    oxybot_lo.plot(lw=2,color='silver',label='LO bottom 20%',ax=axs[2])
    oxytop_orca.plot(lw=1,color='teal',label='ORCA top 20%',ax=axs[2])
    oxybot_orca.plot(lw=1,color='k',label='ORCA bottom 20%',ax=axs[2])
    #axs[2].legend()
    axs[2].grid(True)
    axs[2].set_ylabel('Oxygen [mg/L]')
    axs[2].set_xlabel('Time')

    plt.suptitle(long_name[j])

    fig.set_size_inches(12,10)
    fig.savefig(fn_out+j+'_20.png')