# Vertically integrate the time-averaged vorticity terms.
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        April/2019

import numpy as np
from glob import glob
from os.path import isfile
from netCDF4 import Dataset
from pandas import Timestamp
from xarray import open_dataset, DataArray
from xarray import Dataset as Datasetx


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
headin = '/gpfs/alpine/cli115/proj-shared/apaloczy/vortbdgt_2005-2009avgs/'

# Seasonal and 2005-2009 averages.
fins = ['vortbdgt_2005-2009JFMavg.nc',
        'vortbdgt_2005-2009AMJavg.nc',
        'vortbdgt_2005-2009JASavg.nc',
        'vortbdgt_2005-2009ONDavg.nc',
        'vortbdgt_2005-2009avg.nc']

terms = ['curlnonl', 'betav', 'errcor', 'stretchp', 'errcor', 'curlpgrad', 'curlhdiff', 'curlvdiff', 'res']

fcmsk = '/ccs/home/apaloczy/analysis/data/cache_vortmsks.npz'
fgrd = '/ccs/home/apaloczy/analysis/data/kzit.nc'
dsg = Dataset(fgrd).variables
kmt = dsg['kmt'][:]
fland=kmt==-1 # Land mask.
fcapy = 501
cm2m = 1e-2

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

Terms = dict()
for fin in fins:
    fin = headin + fin
    ds = open_dataset(fin).squeeze(dim='time', drop=True).isel(nlat=slice(fcapy))

    for term in terms:
        vars().update({term:ds[term]})
        print(term)

        # print('Vertically integrating terms.')
        # Vertically integrate terms (interior and topographic budgets).
        Term = vars()[term]
        Iterm = np.nansum(Term*Izimskdzt, axis=0) # [m/s2].
        Iterm[fland] = np.nan

        Tterm = np.nansum(Term*Tzimskdzt, axis=0) # [m/s2].
        Tterm[fland] = np.nan

        Iterm = stripmsk(Iterm)
        Tterm = stripmsk(Tterm)

        iterm = 'I'+term
        tterm = 'T'+term
        coords = dict(lon=(dimsxy, lont), lat=(dimsxy,latt))

        Iterm = DataArray(Iterm, coords=coords,
        dims=dimsxy)
        Tterm = DataArray(Tterm, coords=coords, dims=dimsxy)

        Terms.update({iterm:Iterm})
        Terms.update({tterm:Tterm})

    fout = fin
    Datasetx(data_vars=Terms, coords=coords).to_netcdf(fin.replace('.nc', 'zint.nc'))
