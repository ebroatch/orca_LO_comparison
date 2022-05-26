import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import cmocean.cm as cmo

initials = ['CI', 'PW', 'NB', 'DB', 'HP', 'TW']

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_0 = '/Users/erinbroatch/Documents/Research/Data/ORCA/datasets/'
fn = {}
fn_daily = {}
fn_clim = {}
fn_seas = {}
for j in initials:
    fn[j] = fn_0 + j + '_ds.nc'
    fn_daily[j] = fn_0 + j + '_ds_daily.nc'
    fn_clim[j] = fn_0 + j + '_ds_clim.nc'
    fn_seas[j] = fn_0 + j + '_ds_seas.nc'

ds_daily = {}
ds_clim = {}
ds_seas = {}
for j in initials:
    ds_daily[j] = xr.open_dataset(fn_daily[j])
    ds_clim[j] = xr.open_dataset(fn_clim[j])
    ds_seas[j] = xr.open_dataset(fn_seas[j])

    # Make some sample colormap plots

    #North buoy data availability
    #Make plot
    fig, axs = plt.subplots(3,1,sharex=True)

    #Salinity
    xr.plot.pcolormesh(ds_daily[j]['sal'],'time','depth',cmap=cmo.haline,vmin=22,vmax=32,levels=11,extend='both',ax=axs[0])
    axs[0].invert_yaxis()
    axs[0].set_title('Salinity [psu]')
    axs[0].set_xlabel('')

    #Temperature
    xr.plot.pcolormesh(ds_daily[j]['temp'],'time','depth',cmap=cmo.thermal,vmin=4,vmax=20,levels=9,extend='both',ax=axs[1])
    axs[1].invert_yaxis()
    axs[1].set_title('Temperature [degC]')
    axs[1].set_xlabel('')

    #Oxygen
    xr.plot.pcolormesh(ds_daily[j]['oxy'],'time','depth',cmap=cmo.matter,vmin=0,vmax=20,levels=11,extend='max',ax=axs[2])
    axs[2].invert_yaxis()
    axs[2].set_title('Oxygen [mg/L]')
    axs[2].set_xlabel('Time')

    plt.suptitle(long_name[j])
    plt.show()
    # fig.set_size_inches(12,10)
    # fig.savefig(fn_out+j+'_20.png')

fig, axs = plt.subplots(3,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_daily[j]['sal'],'time','pressure',levels=13,cmap=cmo.haline, vmin=28, vmax=31, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Salinity')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

#Overview of data availability
fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_daily[j]['sal'],'time','pressure',levels=13,cmap=cmo.haline, vmin=28, vmax=31, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Salinity')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_daily[j]['temp'],'time','pressure',levels=16,cmap=cmo.thermal, vmin=5, vmax=20, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Temperature')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_daily[j]['oxy'],'time','pressure',levels=11,cmap=cmo.tempo, vmin=0, vmax=25, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Oxygen')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_daily[j]['par'],'time','pressure',levels=11,cmap=cmo.solar, vmin=0, vmax=50, add_colorbar=False) 
    plt.gca().invert_yaxis()
    plt.ylabel('')
    plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('PAR')
fig.supxlabel('Time')
fig.supylabel('PAR [$\mathrm{\mu Em^{-2} s^{-1}}$]')

fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

#Climatology
fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_clim[j]['sal'],'yearday','pressure',levels=13,cmap=cmo.haline, vmin=28, vmax=31, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Salinity climatology')

fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_clim[j]['temp'],'yearday','pressure',levels=16,cmap=cmo.thermal, vmin=5, vmax=20, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Temperature climatology')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_clim[j]['oxy'],'yearday','pressure',levels=11,cmap=cmo.tempo, vmin=0, vmax=25, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('Oxygen climatology')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

fig, axs = plt.subplots(6,1,sharex=True)
axind = 1
for j in initials:
    plt.subplot(6,1,axind)
    cs=xr.plot.pcolormesh(ds_clim[j]['par'],'yearday','pressure',levels=11,cmap=cmo.solar, vmin=0, vmax=100, add_colorbar=False)
    plt.gca().invert_yaxis()
    if axind !=4:
        plt.ylabel('')
    if axind !=7:
        plt.xlabel('')
    plt.title(long_name[j])
    axind = axind+1
plt.suptitle('PAR climatology')
fig.subplots_adjust(right=0.85, hspace=0.5)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(cs, cax=cbar_ax)

# Seasonal data
fig, axs = plt.subplots(2,3,sharex=True)
axind = 1
for j in initials:
    plt.subplot(2,3,axind)
    plt.plot(ds_seas[j]['sal'].sel(season='DJF'),ds_seas[j]['pressure'],color='tab:blue',label='Winter')
    plt.plot(ds_seas[j]['sal'].sel(season='MAM'),ds_seas[j]['pressure'],color='tab:green',label='Spring')
    plt.plot(ds_seas[j]['sal'].sel(season='JJA'),ds_seas[j]['pressure'],color='tab:red',label='Summer')
    plt.plot(ds_seas[j]['sal'].sel(season='SON'),ds_seas[j]['pressure'],color='tab:orange',label='Autumn')
    plt.title(long_name[j])
    plt.grid(True)
    plt.xlim(24, 32)
    plt.ylim(0,125)
    plt.gca().invert_yaxis()
    if axind==1:
        plt.legend(loc='upper left')
    axind = axind+1
plt.suptitle('Seasonal salinity profiles')
fig.supxlabel('Salinity [psu]')
fig.supylabel('Pressure [dbar]')

fig, axs = plt.subplots(2,3,sharex=True)
axind = 1
for j in initials:
    plt.subplot(2,3,axind)
    plt.plot(ds_seas[j]['temp'].sel(season='DJF'),ds_seas[j]['pressure'],color='tab:blue',label='Winter')
    plt.plot(ds_seas[j]['temp'].sel(season='MAM'),ds_seas[j]['pressure'],color='tab:green',label='Spring')
    plt.plot(ds_seas[j]['temp'].sel(season='JJA'),ds_seas[j]['pressure'],color='tab:red',label='Summer')
    plt.plot(ds_seas[j]['temp'].sel(season='SON'),ds_seas[j]['pressure'],color='tab:orange',label='Autumn')
    plt.title(long_name[j])
    plt.grid(True)
    #plt.xlim(6, 18)
    plt.ylim(0,125)
    plt.gca().invert_yaxis()
    if axind==6:
        plt.legend(loc='lower right')
    axind = axind+1
plt.suptitle('Seasonal temperature profiles')
fig.supxlabel('Temperature [C]')
fig.supylabel('Pressure [dbar]')

fig, axs = plt.subplots(2,3,sharex=True)
axind = 1
for j in initials:
    plt.subplot(2,3,axind)
    plt.plot(ds_seas[j]['oxy'].sel(season='DJF'),ds_seas[j]['pressure'],color='tab:blue',label='Winter')
    plt.plot(ds_seas[j]['oxy'].sel(season='MAM'),ds_seas[j]['pressure'],color='tab:green',label='Spring')
    plt.plot(ds_seas[j]['oxy'].sel(season='JJA'),ds_seas[j]['pressure'],color='tab:red',label='Summer')
    plt.plot(ds_seas[j]['oxy'].sel(season='SON'),ds_seas[j]['pressure'],color='tab:orange',label='Autumn')
    plt.title(long_name[j])
    plt.grid(True)
    plt.xlim(0, 12)
    plt.ylim(0,125)
    plt.gca().invert_yaxis()
    if axind==6:
        plt.legend(loc='lower right')
    axind = axind+1
plt.suptitle('Seasonal oxygen profiles')
fig.supxlabel('Oxygen [mg/L]')
fig.supylabel('Pressure [dbar]')

fig, axs = plt.subplots(2,3,sharex=True)
axind = 1
for j in initials:
    plt.subplot(2,3,axind)
    plt.plot(ds_seas[j]['par'].sel(season='DJF'),ds_seas[j]['pressure'],color='tab:blue',label='Winter')
    plt.plot(ds_seas[j]['par'].sel(season='MAM'),ds_seas[j]['pressure'],color='tab:green',label='Spring')
    plt.plot(ds_seas[j]['par'].sel(season='JJA'),ds_seas[j]['pressure'],color='tab:red',label='Summer')
    plt.plot(ds_seas[j]['par'].sel(season='SON'),ds_seas[j]['pressure'],color='tab:orange',label='Autumn')
    plt.title(long_name[j])
    plt.grid(True)
    plt.xlim(0, 80)
    plt.ylim(0,125)
    plt.gca().invert_yaxis()
    if axind==6:
        plt.legend(loc='lower right')
    axind = axind+1
plt.suptitle('Seasonal PAR profiles')
fig.supxlabel('PAR [$\mathrm{\mu Em^{-2} s^{-1}}$]')
fig.supylabel('Pressure [dbar]')

#Show plots
# plt.show()