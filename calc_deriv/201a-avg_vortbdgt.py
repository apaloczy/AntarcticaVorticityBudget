# Description: Calculate seasonal time averages of the vorticity terms.
#
# Author:      André Palóczy
# E-mail:      paloczy@gmail.com
# Date:        March/2018

import numpy as np
from glob import glob
from os import system
from xarray import open_dataset

def avgds(globpattern, headin):
    fnames = glob(headin+globpattern)
    fnames.sort()
    n=0
    for f in fnames:
        print(f)
        if n>0:
            ds += open_dataset(f)
        else:
            ds = open_dataset(f)
        n+=1

    ds /= n
    return ds

headin = '/gpfs/alpine/cli115/proj-shared/hewang/data/vortbud/SO_curlonly/'
headout = '/gpfs/alpine/cli115/proj-shared/apaloczy/vortbdgt_2005-2009avgs/'
fout_2005_2009 = 'vortbdgt_2005-2009avg.nc'

# Seasonal and 2005-2009 averages.
fouts = ('vortbdgt_2005-2009JFMavg.nc',
         'vortbdgt_2005-2009AMJavg.nc',
         'vortbdgt_2005-2009JASavg.nc',
         'vortbdgt_2005-2009ONDavg.nc')

globpatterns = ('vortbud_????-0[1-3]-??_so_50S.nc',
                'vortbud_????-0[4-6]-??_so_50S.nc',
                'vortbud_????-0[7-9]-??_so_50S.nc',
                'vortbud_????-1[0-2]-??_so_50S.nc')

DS_2005_2009 = None
for globpattern, fout in zip(globpatterns, fouts):
    fnames = glob(headin+globpattern)
    DS = avgds(globpattern, headin)
    DS.to_netcdf(headout+fout)
    if DS_2005_2009 is not None:
        DS_2005_2009 += DS
    else:
        DS_2005_2009 = DS

# Save 2005-2009 average.
DS_2005_2009 /= 4
DS_2005_2009.to_netcdf(headout+fout_2005_2009)
