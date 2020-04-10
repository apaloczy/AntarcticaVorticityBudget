# Description: Calculate daily topographic and interior vorticity terms.
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        March/2018

import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from xarray import open_dataset, DataArray
from xarray import Dataset as Datasetx
from pandas import Timestamp
from os.path import isfile
from glob import glob
# import cmocean.cm as cmo
import cartopy as ctpy
import cartopy.crs as crs
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.ndimage.filters import gaussian_filter


def stripmsk(arr, mask_invalid=False):
    if mask_invalid:
        arr = np.ma.masked_invalid(arr)
    if np.ma.isMA(arr):
        msk = arr.mask
        arr = arr.data
        arr[msk] = np.nan

    return arr


def lon360to180(lon):
	"""
	Converts longitude values in the range [0,360]
	to longitude values in the range [-180,+180].
	"""
	lon = np.asanyarray(lon)
	return ((lon + 180.) % 360.) - 180.




#---
plt.close('all')
head = '/gpfs/alpine/cli115/proj-shared/hewang/data/vortbud/SO_curlonly/'
terms = ['curlnonl', 'betav', 'stretchp', 'errcor', 'curlpgrad', 'curlhdiff', 'curlvdiff', 'res']
# terms = ['betav', 'curlvdiff', 'res']

TOPOG_TERMS = False
PLOT_FIGS = False
SAVE_TERMS_netCDF = True
headout = '/gpfs/alpine/cli115/proj-shared/apaloczy/vortbdgt_2005-2009daily/'

# Get segment lat/lon limits.
segs_lims = {
'Ross':[165., -130., -79., -68.],
'Amundsen-Bellingshausen':[-130., -75., -77., -60.],
'WAP':[-75., -53., -74., -60.],
'Weddell':[-53., -11., -78., -59.],
'W-EA':[-11., 77., -72., -60.],
'E-EA':[77., 165., -72., -60.],
}

fcapy = 501
cm2m = 1e-2
# Load vorticity terms.

# Get isobath coordinates.
bdry_isobs = [800, 2500] # [1000, 2000]
misob = 1000 # Isobath to plot the (T - Tf) section on.
iisob, oisob = np.min(bdry_isobs), np.max(bdry_isobs)

fisobs = '/ccs/home/apaloczy/analysis/data/isobaths.nc'
ncx = Dataset(fisobs)
xi = ncx["%d m isobath (U-points)"%iisob]['xiso'][:]
yi = ncx["%d m isobath (U-points)"%iisob]['yiso'][:]
xm = ncx["%d m isobath (U-points)"%misob]['xiso'][:]
ym = ncx["%d m isobath (U-points)"%misob]['yiso'][:]
xo = ncx["%d m isobath (U-points)"%oisob]['xiso'][:]
yo = ncx["%d m isobath (U-points)"%oisob]['yiso'][:]

fcmsk = '/ccs/home/apaloczy/analysis/data/cache_vortmsks.npz'
# fgrd = '/ccs/home/apaloczy/analysis/data/POP-dzu_dzt_kzit_subsetSO.nc'
fgrd = '/ccs/home/apaloczy/analysis/data/kzit.nc'
fgrd_dz = '/ccs/home/apaloczy/analysis/data/POP-dzu_dzt_kzit_subsetSO.nc'
dsg = Dataset(fgrd).variables
dsg_dzt = Dataset(fgrd_dz).variables
dzt = dsg_dzt['dzt'][:]*cm2m
nz, ny, nx = dzt.shape

# Get z_I and dzu.
kmt = dsg['kmt'][:fcapy,:]
kmtzi = dsg['kzit'][:fcapy,:] # Index of deepest T-cell with no sidewalls.

# Change from Fotran to Python indexing.
kmt = kmt - 1
kmtzi = kmtzi - 1

fland = kmt==-1
kmt[fland] = 0
kmtzi[fland] = 0

if isfile(fcmsk):
    print('Loading cache.')
    d = np.load(fcmsk)
    Izimsk = d['Izimsk']
    Tzimsk = d['Tzimsk']
else:
    print('Calculating Izimsk and Tzimsk.')
    Izimsk = np.zeros((nz, ny, nx)) # Mask for grid points at k > kmtzi.
    Tzimsk = Izimsk.copy()          # Mask for grid points at kmt < k < kmtzi.

    for j in range(ny):
        for i in range(nx):
            ziji = kmtzi[j, i]
            kmtji = kmt[j, i]
            Izimsk[:ziji, j, i] = 1      # Interior mask.
            Tzimsk[ziji:kmtji, j, i] = 1 # Bottom mask.

    Izimsk = np.int32(Izimsk)
    Tzimsk = np.int32(Tzimsk)

    np.savez(fcmsk, Izimsk=Izimsk, Tzimsk=Tzimsk)

lont, latt = dsg_dzt['TLONG'][:], dsg_dzt['TLAT'][:]
lont = lon360to180(lont)
isep = np.abs(np.gradient(lont[0,:])).argmax() + 1

lont = stripmsk(lont)
latt = stripmsk(latt)
dims = ['time', 'y', 'x']
dimsxy = ['y', 'x']

Izimskdzt = Izimsk*dzt
Tzimskdzt = Tzimsk*dzt

IH = Izimskdzt.sum(axis=0)
TH = Tzimskdzt.sum(axis=0)

IH[IH==0] = np.nan
TH[TH==0] = np.nan

fnames = glob(head+'vortbud_20??-??-??_so_50S.nc')
fnames.sort()

# fnames = fnames[-3:]

Terms = dict()
for f in fnames:
    date = f.split('.')[0][-17:-7]
    print(date)

    # print('Loading vorticity terms.')
    ds = open_dataset(f).squeeze(dim='time', drop=True).isel(nlat=slice(fcapy))
    for term in terms:
        vars().update({term:ds[term]})
        # vars().update({term:np.ma.masked_invalid(ds[term].values)})
        # print(term)

        # print('Vertically integrating terms.')
        # Vertically integrate terms (interior and topographic budgets).
        Term = vars()[term]
        # Iterm = np.nansum(Term*Izimskdzt, axis=0)/IH # [1/s2].
        Iterm = np.nansum(Term*Izimskdzt, axis=0) # [m/s2].
        Iterm[fland] = np.nan

        if TOPOG_TERMS:
            # Tterm = np.nansum(Term*Tzimskdzt, axis=0)/TH # [1/s2].
            Tterm = np.nansum(Term*Tzimskdzt, axis=0) # [m/s2].
            Tterm[fland] = np.nan

        Iterm = stripmsk(Iterm)
        if TOPOG_TERMS:
            Tterm = stripmsk(Tterm)

        iterm = 'I'+term
        if TOPOG_TERMS:
            tterm = 'T'+term
        # Add timestamp.

        coords = dict(lon=(dimsxy, lont), lat=(dimsxy,latt))

        Iterm = DataArray(Iterm, coords=coords,
        dims=dimsxy)
        if TOPOG_TERMS:
            Tterm = DataArray(Tterm, coords=coords, dims=dimsxy)

        t = np.array(Timestamp(date))
        Iterm.coords.update(dict(time=t))

        Terms.update({iterm:Iterm})
        if TOPOG_TERMS:
            Terms.update({tterm:Tterm})

    if SAVE_TERMS_netCDF:
        fout = headout + 'vortbdgt_' + date + '.nc'
        Datasetx(data_vars=Terms, coords=coords).to_netcdf(fout, unlimited_dims='time')
