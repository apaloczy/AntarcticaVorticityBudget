# Description: Calculate total ocean stress curl (ice + wind).
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        April/2019

import numpy as np
# import matplotlib.pyplot as plt
from glob import glob
from pandas import Timestamp
from xarray import open_dataset, DataArray, Dataset


def zcurl(Ux, Uy, dxu, dyu, tarea, kmt, cyclic=True):
    """
    USAGE
    -----
    curlzU = zcurl(Ux, Uy, dxu, dyu, tarea, kmt, cyclic=True)
    z-curl operator on vector fields defined at U-points:

    curlzU(Ux,Uy) =  (1/dy)*Dx[Ay(dy*Uy)] - (1/dx)*Dy[Ax(dx*Ux)]

    This routine returns curl at T-points.

    Based on POP subroutine "zcurl" in operators.F90
    Author: Steve Yeager (NCAR).
    Translated from NCL to Python by André Palóczy (SIO/UCSD).
    """
    p5 = 0.5
    assert Ux.shape==Uy.shape, "Ux and Uy must have the same shape."
    ndims = Ux.ndim
    assert Ux.ndim==2 and Uy.ndim==2, "Inputs must be (x,y) 2D fields."
    ny, nx = Ux.shape

    worky = p5*Uy*dyu
    workx = p5*Ux*dxu
    curlzU = workx.copy()

    workxNE = workx.copy()                   # [j,i] term
    workxSE = workx.copy()
    workxSE[1:-1, :] = workx[:-2, :]         # [j-1, i] term
    workxSE[0, :] = np.nan
    workxNW = workx.copy()
    workxNW[:, 1:-1] = workx[:, :-2]         # [j, i-1] term
    if cyclic:
        workxNW[:, 0] = workx[:, -1]         # wraparound
    else:
        workxNW[:, 0] = np.nan
    workxSW = workx.copy()
    workxSW[1:-1, 1:-1] = workx[:-2, :-2]    # [j-1, i-1] term
    if cyclic:
        workxSW[1:-1, 0] = workx[:-2, -1]    # wraparound
    else:
        workxSW[1:-1, 0] = np.nan
    workxSW[0, :] = np.nan

    workyNE = worky.copy()                   # [j,i] term
    workySE = worky.copy()
    workySE[1:-1, :] = worky[:-2, :]         # [j-1, i] term
    workySE[0, :] = np.nan
    workyNW = worky.copy()
    workyNW[:, 1:-1] = worky[:, :-2]         # [j, i-1] term
    if cyclic:
        workyNW[:, 0] = worky[:, -1]         # wraparound
    else:
        workyNW[:, 0] = np.nan
    workySW = worky.copy()
    workySW[1:-1, 1:-1] = worky[:-2, :-2]    # [j-1, i-1] term
    if cyclic:
        workySW[1:-1, 0] = worky[:-2, -1]    # wraparound
    else:
        workySW[1:-1, 0] = np.nan
    workySW[0, :] = np.nan

    # d(Uy)/dx, -d(Ux)/dy
    curlzU = ((workyNE+workySE)-(workyNW+workySW)) - ((workxNE+workxNW)-(workxSE+workxSW))

    curlzU = curlzU/tarea  # Divide by area.
    curlzU[0, :] = np.nan  # Curl should be undefined at j=0 T-points.

    return curlzU


#---
# plt.close('all')

fname_out = 'srfstresses_ocn.nc'
fname_out_monthly = 'srfstresses_ocn_monthly.nc'
fname_out_ssn = 'srfstresses_ocn_seasonal.nc'

head_out = '/gpfs/alpine/cli115/proj-shared/apaloczy/wnd-ice-ocn_stress/'
YEAR_START = 2005#1959#1969
cgs2mks_stress = 0.1 # dyne/cm2 to N/m2.
cm2m = 1e-2 # cm to m.

head_fin = '/gpfs/alpine/cli115/proj-shared/ia_top_tx0.1_v2_60yrs/'
tail_file_ocn = 'ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.????-??.nc'
tail_file_ice = 'ia_top_tx0.1_v2_yel_patc_1948_intel.cice.h.????-??.nc'
fgrd = '~/analysis/data/POP-dzu_dzt_kzit_subsetSO.nc'
dsg = open_dataset(fgrd)
kmt = dsg['kmt'][:]
fland=kmt==-1 # Land mask.

fdirs_ocn = glob(head_fin+'ia_top_tx0.1_v2_yel_patc_1948_intel_def_year_????/ocn/hist/')
fdirs_ocn.sort()
fdirs_ice = glob(head_fin+'ia_top_tx0.1_v2_yel_patc_1948_intel_def_year_????/ice/hist/')
fdirs_ice.sort()
if not isinstance(fdirs_ocn, list):
    fdirs_ocn = [fdirs_ocn]
    fdirs_ice = [fdirs_ice]
fnames_ocn = []
for fdir in fdirs_ocn:
    if int(fdir.split('/ocn/')[0][-4:])<YEAR_START:
        continue
    fnamesi = glob(fdir + tail_file_ocn)
    fnamesi.sort()
    for f in fnamesi:
        fnames_ocn.append(f)
fnames_ice = []
for fdir in fdirs_ice:
    if int(fdir.split('/ice/')[0][-4:])<YEAR_START:
        continue
    fnamesi = glob(fdir + tail_file_ice)
    fnamesi.sort()
    for f in fnamesi:
        fnames_ice.append(f)

nt_ocn = len(fnames_ocn) # Total number of files.
nt_ice = len(fnames_ice)
assert nt_ocn==nt_ice
nt = nt_ocn

fcap = 501

for n in range(nt):
    print("File ",n+1,"/",nt)
    ocn = open_dataset(fnames_ocn[n])
    ice = open_dataset(fnames_ice[n])
    if n==0:
        lonu = ocn['ULONG'][:fcap,:].values
        latu = ocn['ULAT'][:fcap,:].values
        dxu = ocn['DXU'][:fcap,:].values*cm2m
        dyu = ocn['DYU'][:fcap,:].values*cm2m
        tarea = ocn['TAREA'][:fcap,:].values*cm2m*cm2m

    taux_wnd = np.squeeze(ocn['TAUX'])[:fcap,:].values*cgs2mks_stress
    tauy_wnd = np.squeeze(ocn['TAUY'])[:fcap,:].values*cgs2mks_stress
    taux_ice = np.squeeze(ice['strocnx'])[:fcap,:].values
    tauy_ice = np.squeeze(ice['strocny'])[:fcap,:].values

    # Calculate ice concentration and total ocean stress.
    # taux,y_ice has already been weighted by the sea ice concentration in the cell.
    iceconc = np.squeeze(ice['aice'])[:fcap,:].values/100
    # taux_tot = iceconc*taux_ice + (1-iceconc)*taux_wnd
    # tauy_tot = iceconc*tauy_ice + (1-iceconc)*tauy_wnd

    taux_tot = taux_wnd - taux_ice
    tauy_tot = tauy_wnd - tauy_ice

    # taux_tot = taux_wnd + taux_ice
    # tauy_tot = tauy_wnd + tauy_ice

    # Want -(ice stress) because it is the stress applied at the surface of the OCEAN.
    # The quantity saved in the CICE outputs is the ocean stress ON THE SEA ICE.

    # Calculate stress curls.
    # curltau_wnd = zcurl(taux_wnd, tauy_wnd, dxu, dyu, tarea, kmt)
    # curltau_ice = zcurl(taux_ice, tauy_ice, dxu, dyu, tarea, kmt)
    curltau_tot = zcurl(taux_tot, tauy_tot, dxu, dyu, tarea, kmt)

    if n==0:
        Taux_tot = taux_tot
        Tauy_tot = tauy_tot
        # Taux_wnd = taux_wnd
        # Tauy_wnd = tauy_wnd
        # Taux_ice = taux_ice
        # Tauy_ice = tauy_ice
        # Curltau_wnd = curltau_wnd
        # Curltau_ice = curltau_ice
        Curltau_tot = curltau_tot
        Iceconc = iceconc
    else:
        Taux_tot = np.dstack((Taux_tot, taux_tot[...,np.newaxis]))
        Tauy_tot = np.dstack((Tauy_tot, tauy_tot[...,np.newaxis]))
        # Taux_wnd = np.dstack((Taux_wnd, taux_wnd[...,np.newaxis]))
        # Tauy_wnd = np.dstack((Tauy_wnd, tauy_wnd[...,np.newaxis]))
        # Taux_ice = np.dstack((Taux_ice, taux_ice[...,np.newaxis]))
        # Tauy_ice = np.dstack((Tauy_ice, tauy_ice[...,np.newaxis]))
        # Curltau_wnd = np.dstack((Curltau_wnd, curltau_wnd[...,np.newaxis]))
        # Curltau_ice = np.dstack((Curltau_ice, curltau_ice[...,np.newaxis]))
        Curltau_tot = np.dstack((Curltau_tot, curltau_tot[...,np.newaxis]))
        Iceconc = np.dstack((Iceconc, iceconc[...,np.newaxis]))

# Make time axis.
t = np.array([])
for fname in fnames_ocn:
    date = fname.split('.')[-2]
    t = np.append(t, Timestamp(date+" 15 12:00:00").to_pydatetime())

dims = ('x', 'y', 't')
coords = dict(lon=(['x', 'y'], lonu), lat=(['x', 'y'], latu), t=('t', t))
Taux_tot = DataArray(Taux_tot, coords=coords, dims=dims)
Tauy_tot = DataArray(Tauy_tot, coords=coords, dims=dims)
# Taux_wnd = DataArray(Taux_wnd, coords=coords, dims=dims)
# Tauy_wnd = DataArray(Tauy_wnd, coords=coords, dims=dims)
# Taux_ice = DataArray(Taux_ice, coords=coords, dims=dims)
# Tauy_ice = DataArray(Tauy_ice, coords=coords, dims=dims)
# Curltau_wnd = DataArray(Curltau_wnd, coords=coords, dims=dims)
# Curltau_ice = DataArray(Curltau_ice, coords=coords, dims=dims)
Curltau_tot = DataArray(Curltau_tot, coords=coords, dims=dims)
Iceconc = DataArray(Iceconc, coords=coords, dims=dims)

# data_vars = dict(tauxwnd=Taux_wnd, tauywnd=Tauy_wnd, tauxice=Taux_ice, tauyice=Tauy_ice, curlwind=Curltau_wnd, curlice=Curltau_ice, curltot=Curltau_tot)
data_vars = dict(tauxtot=Taux_tot, tauytot=Tauy_tot, curltot=Curltau_tot, iceconc=Iceconc)
ds = Dataset(data_vars=data_vars, coords=coords)

# Mothly fields.
ds.to_netcdf(head_out+fname_out)
ds.close()

# Monthly averages.
dsm = ds.groupby('t.month').mean(dim='t')
dsm.to_netcdf(head_out+fname_out_monthly)
dsm.close()

# Seasonal averages.
dsssn = ds.groupby('t.season').mean(dim='t')
dsssn.to_netcdf(head_out+fname_out_ssn)
dsssn.close()
