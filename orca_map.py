import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import os
import sys
import cmocean.cm as cmo
import datetime as dt
import netCDF4 as nc

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

# Get useful LiveOcean functions
pth=os.path.abspath('/Users/erinbroatch/Documents/GitHub/LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zrfun

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

def add_coast(ax, dir0, color='k'):
    fn = dir0 + 'coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    ax.plot(C['lon'].values, C['lat'].values, '-', color=color, linewidth=0.5)

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
im2=ax.pcolormesh(lon,lat,waterplot,cmap=cmo.ice_r,vmin=0, vmax=300, shading='auto')
ax.set_xlim(-123.3, -122)
ax.set_ylim(46.9, 48.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
dar(ax)
add_coast(ax,dir0='/Users/erinbroatch/Documents/Research/Data/')
ax.set_title('Puget Sound')
#cbar = fig.colorbar(im2, ax=ax, label='Depth $\mathrm{[m]}$')

# Add moorings
for j in initials:
    ax.scatter(lonlat[j][0],lonlat[j][1],c=plot_colors[j],s=25,label=long_name[j])
    ax.text(lonlat[j][0]+0.02,lonlat[j][1],j,c=plot_colors[j],fontsize=13, fontweight='bold')
plt.legend(loc='lower right')

# Tidy
ax.set_title('ORCA mooring locations', fontsize=13, fontweight=5)
fig.set_size_inches(8,10)

# Save plot
plt.savefig(fn_out+'orca_map.png', dpi=300)

