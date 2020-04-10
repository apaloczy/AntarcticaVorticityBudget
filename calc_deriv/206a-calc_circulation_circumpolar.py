# Description: Calculate coherences between daily vorticity terms.
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        April/2019

import numpy as np
from glob import glob
from netCDF4 import Dataset
from xarray import open_dataset, DataArray
from xarray import Dataset as Datasetx
from pandas import Timestamp


def lon360to180(lon):
	"""
	Converts longitude values in the range [0,360]
	to longitude values in the range [-180,+180].
	"""
	lon = np.asanyarray(lon)
	return ((lon + 180.) % 360.) - 180.


#---
head = '/ccs/home/apaloczy/analysis/'
headin = '/gpfs/alpine/cli115/proj-shared/apaloczy/vortbdgt_2005-2009daily/'
terms = ['Ibetav', 'Icurlvdiff', 'Icurlhdiff', 'Istretchp', 'Ires', 'Icurlnonl']

fnames = glob(headin+'*.nc')
fnames.sort()

# Get isobath masks.
fouter = head+'data/volmsk2500m.npz'
finner = head+'data/volmsk800m.npz'
fname_ncgrid = head+'data/POP-dzu_dzt_kzit_subsetSO.nc'
cm2m = 1e-2

ncg = Dataset(fname_ncgrid)
latt = ncg.variables['TLAT'][:]
lont = lon360to180(ncg.variables['TLONG'][:])
TAREA = Dataset(head+'data/TAREA.nc')['TAREA'][:501,:]*cm2m**2

d = np.load(fouter)
xo, yo, msko = d['lont'], d['latt'], np.bool8(d['volmsk'])
d = np.load(finner)
xi, yi, mski = d['lont'], d['latt'], np.bool8(d['volmsk'])
msk = np.logical_and(msko, ~mski)
TAREAmsk = TAREA[msk]
A = TAREAmsk.sum()

# fnames = fnames[-90:]
t = np.array([])
Cterms = dict()
_ = [Cterms.update({term:np.array([], ndmin=1)}) for term in terms]
n=0
for fname in fnames:
	date = fname.split('.')[0][-10:]
	print(date)
	t = np.append(t, Timestamp(date).to_pydatetime())
	ds = open_dataset(fname)
	for term in terms:
		dsterm = ds[term]
		dsterm[:,0] = np.nan
		dsterm[:,-1] = np.nan
		aux = np.nansum(dsterm.values[msk]*TAREAmsk)/A # Area-averaging.
		# NOTE: Excluding first and last columns because curl did not wrap around.
		Cterms[term] = np.append(Cterms[term], aux)

for term in terms:
	Cterms[term] = DataArray(Cterms[term], coords=dict(t=t), dims='t')

ds = Datasetx(data_vars=Cterms, coords=dict(t=t))

fout = head+'data/circulation_terms_circumpolar.nc'
ds.to_netcdf(fout)
