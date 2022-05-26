import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import cmocean.cm as cmo
import datetime as dt
import netCDF4 as nc
from lo_tools import zrfun #LOENV MUST BE ACTIVATED
from lo_tools import zfun
from lo_tools import plotting_functions as plotfun

initials = ['CI', 'PW', 'NB', 'DB', 'HP', 'TW']

plot_colors = {'CI':'tab:blue',
    'PW':'tab:pink',
    'NB':'tab:green',
    'DB':'tab:red',
    'HP':'tab:orange',
    'TW':'tab:purple'}

long_name = {'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'North Buoy (Hansville)',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'}

fn_out = '/Users/erinbroatch/Documents/orca_report/figures/'

lonlat = {'CI': [-122.7300, 47.2800],
    'PW': [-122.3972, 47.7612],
    'NB': [-122.6270, 47.9073],
    'DB': [-122.8029, 47.8034],
    'HP': [-123.1126, 47.4218],
    'TW': [-123.0083, 47.3750]}

# Set plot style
plt.close('all')
plt.style.use('seaborn-darkgrid')

# Get data for bathy map
fn='/Users/erinbroatch/Documents/Research/Output/2017_01_19_1900.nc'
ds = nc.Dataset(fn)
zeta = ds['zeta'][0,:,:]
ds.close()
# Get grid info
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = G['h']
hplot=np.ma.masked_array(h,mask=zeta.mask)
lon = G['lon_rho']
lat = G['lat_rho']
waterplot=np.ma.masked_array(np.zeros(h.shape),mask=zeta.mask)

############
# Make plot
fig, ax = plt.subplots(1,1)

# Make bathy and coastline plot
im2=ax.pcolormesh(lon,lat,hplot,cmap=cmo.ice_r,vmin=0, vmax=200, shading='auto')
ax.set_xlim(-122.725, -122.55)
ax.set_ylim(47.85, 47.95)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plotfun.dar(ax)
plotfun.add_coast(ax,dir0='/Users/erinbroatch/Documents/Research/Data/')
ax.set_title('Puget Sound')
cbar = fig.colorbar(im2, ax=ax, label='Depth $\mathrm{[m]}$')

# Add mooring
ax.scatter(lonlat['NB'][0],lonlat['NB'][1],c='tab:red',s=25,label=long_name['NB'])
ax.text(lonlat['NB'][0]+0.002,lonlat['NB'][1],'NB',c='tab:red',fontsize=13, fontweight='bold')

# Tidy
ax.set_title('North Buoy adcp', fontsize=13, fontweight=5)
fig.set_size_inches(10,8)

# Save plot
plt.savefig(fn_out+'adcp_map.png', dpi=300)