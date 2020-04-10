# Calculate KZIT (index of first interior grid cell)
# interior grid cell = cell with no adjacent sidewall boundaries.
#
# Minimum of the 8 adjacent T-cells and the current T cell, -1 level.
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        May/2019

import numpy as np
from netCDF4 import Dataset
from xarray import DataArray
from xarray import Dataset as Datasetx

##
continent_flag = -1
fname = "/lustre/atlas1/cli115/proj-shared/ia_top_tx0.1_v2_60yrs/ia_top_tx0.1_v2_yel_patc_1948_intel_def_year_2009/ocn/hist/ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.2009-01.nc"

nc = Dataset(fname)

lont = nc.variables['TLONG'][:]
latt = nc.variables['TLAT'][:]
z_t = nc.variables['z_t'][:]
h_t = nc.variables['HT'][:]
kmt = nc.variables['KMT'][:] # k index of deepest T-cell.

ny, nx = h_t.shape
nym = ny - 1
nxm = nx - 1
nz = z_t.size
kmt = kmt - 1 # Python indexing starts at 0, FORTRAN indexing starts at 1.

# Get index of deepest interior grid cell.
print("Calculating KZIT.")
kzit = np.zeros((ny,nx))

for j in range(ny):
    print(j+1," / ",ny)
    for i in range(nx):
        if kmt[j,i]==-1: # Continent mask.
            kzit[j,i] = -1
            continue
        else:
            if np.logical_and(j==0, i==0):
                kzit[j,i] = np.min([kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,i+1], kmt[j,i+1], kmt[j,i]]) - 1
            elif np.logical_and(j==nym, i==0):
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j,i+1], kmt[j-1,i+1], kmt[j-1,i], kmt[j,i]]) - 1
            elif np.logical_and(j==nym, i==nxm):
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j,0], kmt[j-1,0], kmt[j-1,i], kmt[j,i]]) - 1
            elif np.logical_and(j==0, i==nxm):
                kzit[j,i] = np.min([kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,0], kmt[j,0], kmt[j,i]]) - 1
            elif j==0:
                kzit[j,i] = np.min([kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,i+1], kmt[j,i+1], kmt[j,i]]) - 1
            elif j==nym:
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j,i+1], kmt[j-1,i+1], kmt[j-1,i], kmt[j,i]]) - 1
            elif i==0:
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,i+1], kmt[j,i+1], kmt[j-1,i+1], kmt[j-1,i], kmt[j,i]]) - 1
            elif i==nxm:
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,0], kmt[j,0], kmt[j-1,0], kmt[j-1,i], kmt[j,i]]) - 1
            else:
                kzit[j,i] = np.min([kmt[j-1,i-1], kmt[j,i-1], kmt[j+1,i-1], kmt[j+1,i], kmt[j+1,i+1], kmt[j,i+1], kmt[j-1,i+1], kmt[j-1,i], kmt[j,i]]) - 1


# Change back to fortran indexing. Continent mask is now 0.
kzit = np.int32(kzit + 1)
kmt = np.int32(kmt + 1)

kzit[kzit==-1] = 0

dims = ['x', 'y']
coords = dict(lont=(dims, lont), latt=(dims, latt))
kmt = DataArray(kmt, coords=coords, dims=dims)
kzit = DataArray(kzit, coords=coords, dims=dims)

fout = '/ccs/home/apaloczy/analysis/data/kzit.nc'
Datasetx(data_vars=dict(kmt=kmt, kzit=kzit), coords=coords).to_netcdf(fout)
