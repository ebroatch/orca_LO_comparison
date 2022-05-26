import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import datetime as dt
import netCDF4 as nc

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
fn_out = '/Users/erinbroatch/Documents/orca_report/figures/'

for j in initials:

    #Salinity plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_daily[j]['sal'],'time','depth',cmap=cmo.haline,vmin=22,vmax=32,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['sal'],'time','depth',cmap=cmo.haline,vmin=22,vmax=32,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['sal']-ds_daily[j]['sal'],'time','depth',cmap=cmo.balance,center=0,vmin=-5,vmax=5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$S [psu]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2017,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' Salinity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_sal.png')

    #Temperature plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_daily[j]['temp'],'time','depth',cmap=cmo.thermal,vmin=4,vmax=20,levels=9,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['temp'],'time','depth',cmap=cmo.thermal,vmin=4,vmax=20,levels=9,extend='both')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['temp']-ds_daily[j]['temp'],'time','depth',cmap=cmo.balance,center=0,vmin=-5,vmax=5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$T [degC]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2017,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' Temperature')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_temp.png')

    #Oxygen plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_daily[j]['oxy'],'time','depth',cmap=cmo.matter,vmin=0,vmax=20,levels=11,extend='max')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['oxy'],'time','depth',cmap=cmo.matter,vmin=0,vmax=20,levels=11,extend='max')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['oxy']-ds_daily[j]['oxy'],'time','depth',cmap=cmo.balance,center=0,vmin=-10,vmax=10,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$ Oxygen [mg/L]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2017,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' Oxygen')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_oxy.png')


