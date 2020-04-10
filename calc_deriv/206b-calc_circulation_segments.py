# Description: Calculate coherences between daily vorticity terms.
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        April/2019

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from netCDF4 import Dataset
from xarray import open_dataset
from pandas import Timestamp
from xarray import DataArray
from xarray import Dataset as Datasetx

def lon360to180(lon):
	"""
	Converts longitude values in the range [0,360]
	to longitude values in the range [-180,+180].
	"""
	lon = np.asanyarray(lon)
	return ((lon + 180.) % 360.) - 180.


# Get segment lat/lon limits.
segs_lims = {
'Ross':[165., -130., -79., -68.],
'Amundsen-Bellingshausen':[-130., -75., -75., -68.],
'WAP':[-75., -53., -74., -60.],
'Weddell':[-62., -11., -79., -59.],
'W-EA':[-11., 77., -72., -60.],
'E-EA':[77., 165., -72., -60.],
}
segments = segs_lims.keys()

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


Segments = dict()
for segment in segments:
	bbox = segs_lims[segment]
	if segment=='Ross':
		in_seg = np.logical_or(lont>=bbox[0], lont<bbox[1])
	else:
		in_seg = np.logical_and(lont>=bbox[0], lont<bbox[1])
	mskn = np.logical_and(msk, in_seg)
	TAREAn = TAREA[mskn]
	An = TAREAn.sum()
	# print(An)
	# plt.figure()
	# plt.imshow(mskn)

	t = np.array([])
	Cterms = dict()
	_ = [Cterms.update({term:np.array([], ndmin=1)}) for term in terms]
	for fname in fnames:
		date = fname.split('.')[0][-10:]
		print(date)
		t = np.append(t, Timestamp(date).to_pydatetime())
		ds = open_dataset(fname)
		for term in terms:
			print(np.isnan(ds[term].values[mskn]).sum())
			dsterm = ds[term]
			dsterm[:,0] = np.nan
			dsterm[:,-1] = np.nan
			aux = np.nansum(dsterm.values[mskn]*TAREAn)/An # Area-averaging.
			Cterms[term] = np.append(Cterms[term], aux)
	Segments.update({segment:Cterms})

for segment in segments:
	Cterms = Segments[segment]
	for term in terms:
		Cterms[term] = DataArray(Cterms[term], coords=dict(t=t), dims='t')

	fout = head+'data/circulation_terms-%s.nc'%segment
	ds = Datasetx(data_vars=Cterms, coords=dict(t=t))
	ds.to_netcdf(fout)
