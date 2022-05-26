import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import datetime as dt
import netCDF4 as nc

import os
import sys
pth=os.path.abspath('/Users/erinbroatch/Documents/GitHub/LO/lo_tools/lo_tools')
if pth not in sys.path:
    sys.path.append(pth)
import zfun

def rot_vec(u,v,theta):
    ur = u*np.cos(theta) + v*np.sin(theta)
    vr = v*np.cos(theta) - u*np.sin(theta)
    return ur, vr

plt.close('all')

initials = ['NB']

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

# ORCA adcp data
fn_0 = '/Users/erinbroatch/Documents/Research/Data/orca_adcp/datasets/'
fn = {}
for j in initials:
    fn[j] = fn_0 + j + '_adcp_ds.nc'
ds = {}
for j in initials:
    ds[j] = xr.open_dataset(fn[j])

fn_hourly = {}
for j in initials:
    fn_hourly[j] = fn_0 + j + '_adcp_ds_hourly.nc'
ds_hourly = {}
for j in initials:
    ds_hourly[j] = xr.open_dataset(fn_hourly[j])

# LO velocity data
fn_lo_0 = '/Users/erinbroatch/Documents/Research/Data/LO_orca_moor/datasets/'
fn_lo = {}
#can add more versions here
for j in initials:
    fn_lo[j] = fn_lo_0 + j + '_adcp_ds.nc'
ds_lo = {}
for j in initials:
    ds_lo[j] = xr.open_dataset(fn_lo[j])

# Location to save figures
fn_out = '/Users/erinbroatch/Documents/Research/Figures/2022_February_8/adcp/'

for j in initials:

    #u plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_hourly[j]['u'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['u'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['u']-ds_hourly[j]['u'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$u [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' u velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_u.png')

    #v plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_hourly[j]['v'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['v'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['v']-ds_hourly[j]['v'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$v [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' v velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_v.png')

    #w plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(ds_hourly[j]['w'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.05,vmax=0.05,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(ds_lo[j]['w'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.05,vmax=0.05,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(ds_lo[j]['w']-ds_hourly[j]['w'],'time','depth',cmap=cmo.balance,center=0,vmin=-0.05,vmax=0.05,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$w [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO-ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' w velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_w.png')

    #echo plot
    fig, axs = plt.subplots(1,1)
    plt.subplot(1,1,1)
    cs=xr.plot.pcolormesh(ds_hourly[j]['echo'],'time','depth',cmap=cmo.dense)#,center=0,vmin=-0.05,vmax=0.05,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' echo')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,5)
    fig.savefig(fn_out+j+'_echo.png')

    #Rotate
    ubar_orca = ds_hourly[j]['u'].mean(dim='z')
    vbar_orca = ds_hourly[j]['v'].mean(dim='z')
    up_orca = ubar_orca - ubar_orca.mean()
    vp_orca = vbar_orca - vbar_orca.mean()
    theta_orca = 0.5 * np.arctan2(2*np.nanmean(up_orca*vp_orca),(np.nanvar(up_orca)-np.nanvar(vp_orca))) #pca, might change
    ubar_orca_rot, vbar_orca_rot = rot_vec(ubar_orca, vbar_orca, theta_orca)
    u_orca_rot, v_orca_rot = rot_vec(ds_hourly[j]['u'],ds_hourly[j]['v'],theta_orca)

    ubar_lo = ds_lo[j]['u'].sel(time=slice('2021-01-01', '2021-12-31')).mean(dim='z')
    vbar_lo = ds_lo[j]['v'].sel(time=slice('2021-01-01', '2021-12-31')).mean(dim='z')
    up_lo = ubar_lo - ubar_lo.mean()
    vp_lo = vbar_lo - vbar_lo.mean()
    theta_lo = 0.5 * np.arctan2(2*np.nanmean(up_lo*vp_lo),(np.nanvar(up_lo)-np.nanvar(vp_lo))) #pca, might change
    ubar_lo_rot, vbar_lo_rot = rot_vec(ubar_lo, vbar_lo, theta_lo)
    u_lo_rot, v_lo_rot = rot_vec(ds_lo[j]['u'].sel(time=slice('2021-01-01', '2021-12-31')),ds_lo[j]['v'].sel(time=slice('2021-01-01', '2021-12-31')),theta_lo)

    #Filter
    u_orca_filt = zfun.lowpass(np.transpose(u_orca_rot.values), f='godin')
    v_orca_filt = zfun.lowpass(np.transpose(v_orca_rot.values), f='godin')
    u_lo_filt = zfun.lowpass((u_lo_rot.values), f='godin') #NEED TO FIX ORDERING
    v_lo_filt = zfun.lowpass((v_lo_rot.values), f='godin')

    u_orca_filt_da = xr.DataArray(u_orca_filt, dims=('time', 'z'), coords={'time':('time',u_orca_rot.time.data), 'depth':('z', u_orca_rot.depth.data)})
    v_orca_filt_da = xr.DataArray(v_orca_filt, dims=('time', 'z'), coords={'time':('time',v_orca_rot.time.data), 'depth':('z', v_orca_rot.depth.data)})
    u_lo_filt_da = xr.DataArray(u_lo_filt, dims=('time', 'z'), coords={'time':('time',u_lo_rot.time.data), 'depth':('z', u_lo_rot.depth.data)})
    v_lo_filt_da = xr.DataArray(v_lo_filt, dims=('time', 'z'), coords={'time':('time',v_lo_rot.time.data), 'depth':('z', v_lo_rot.depth.data)})

    #Ellipse scatter
    fig, axs = plt.subplots(1,2,sharey=True)
    plt.subplot(1,2,1)
    plt.scatter(ubar_orca,vbar_orca,s=50,alpha=0.2,edgecolors='None')
    plt.title('ORCA')
    plt.xlabel('u')
    plt.ylabel('v')
    plt.xlim(-0.4,0.4)
    plt.ylim(-0.4,0.4)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    plt.subplot(1,2,2)
    plt.scatter(ubar_lo,vbar_lo,s=50,alpha=0.2,edgecolors='None')
    plt.title('LO')
    plt.xlabel('u')
    plt.xlim(-0.4,0.4)
    plt.ylim(-0.4,0.4)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    fig.set_size_inches(10,5)
    fig.savefig(fn_out+j+'_scatter.png')

    # Rotated ellipse scatter
    fig, axs = plt.subplots(1,2,sharey=True)
    plt.subplot(1,2,1)
    plt.scatter(ubar_orca_rot,vbar_orca_rot,s=50,alpha=0.2,edgecolors='None')
    plt.title('ORCA, theta='+str(theta_orca))
    plt.xlabel('u')
    plt.ylabel('v')
    plt.xlim(-0.6,0.6)
    plt.ylim(-0.6,0.6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    plt.subplot(1,2,2)
    plt.scatter(ubar_lo_rot,vbar_lo_rot,s=50,alpha=0.2,edgecolors='None')
    plt.title('LO, theta='+str(theta_lo))
    plt.xlabel('u')
    plt.xlim(-0.6,0.6)
    plt.ylim(-0.6,0.6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    fig.set_size_inches(10,5)
    fig.savefig(fn_out+j+'_scatter_rot.png')

#u rot plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(u_orca_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(-u_lo_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('-LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(u_lo_rot+u_orca_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$u [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO+ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' along-channel velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_u_rot.png')

    #v rot plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(v_orca_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(-v_lo_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('-LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(v_lo_rot+v_orca_rot,'time','depth',cmap=cmo.balance,center=0,vmin=-0.5,vmax=0.5,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$v [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO+ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' cross-channel velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_v_rot.png')

#u filt plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(u_orca_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.3,vmax=0.3,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(-u_lo_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.3,vmax=0.3,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('-LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(u_lo_filt_da+u_orca_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.3,vmax=0.3,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$u [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO+ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+'filtered along-channel velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_u_filt.png')

    #v filt plot
    fig, axs = plt.subplots(3,1,sharex=True)
    plt.subplot(3,1,1)
    cs=xr.plot.pcolormesh(v_orca_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.2,vmax=0.2,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('ORCA')
    plt.xlabel('')

    plt.subplot(3,1,2)
    cs=xr.plot.pcolormesh(-v_lo_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.2,vmax=0.2,levels=11,extend='both')
    plt.gca().invert_yaxis()
    plt.title('-LO')
    plt.xlabel('')

    plt.subplot(3,1,3)
    cs=xr.plot.pcolormesh(v_lo_filt_da+v_orca_filt_da,'time','depth',cmap=cmo.balance,center=0,vmin=-0.2,vmax=0.2,levels=11,extend='both',cbar_kwargs={'label':'$\Delta$v [m/s]'})
    plt.gca().invert_yaxis()
    plt.title('LO+ORCA')

    plt.xlim([dt.datetime(2021,1,1),dt.datetime(2021,12,31)])
    plt.suptitle(long_name[j]+' filtered cross-channel velocity')
    fig.subplots_adjust(right=0.85, hspace=0.5)

    fig.set_size_inches(10,10)
    fig.savefig(fn_out+j+'_v_filt.png')
