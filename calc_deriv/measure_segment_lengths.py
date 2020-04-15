# Measure length of each segment (except WAP)
# from longitude limits.
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def near(x, x0):
    return np.abs(x - x0).argmin()

#---
plt.close('all')

head_data = "../data_reproduce_figs/"
isob = 800

# Get segment lat/lon limits.
segs_lims = {
'Amundsen-Bellingshausen':[-130., -75., -75., -68.],
'WAP':[-75., -55., -74., -60.],
'Weddell':[-62., -11., -79., -59.],
'Ross':[345., 410., -79., -68.],
'W-EA':[-11., 77., -71., -60.],
'E-EA':[77., 165., -71., -60.],
}

ds = Dataset(head_data+"isobaths.nc")['%d m isobath (U-points)'%isob]
dist = ds.variables['diso'][:]
long = ds.variables['xiso'][:]

fl, fr = near(long, -130), near(long, -75)
Length_AB = dist[fr] - dist[fl]

fl, fr = near(long, -75), near(long, -55)
Length_WAP = dist[fr] - dist[fl]

fl, fr = near(long, -62), near(long, -11)
Length_Weddell = dist[fr] - dist[fl]

fl, fr = near(long, -11), near(long, 77)
Length_WEA = dist[fr] - dist[fl]

fl, fr = near(long, 77), near(long, 165)
Length_EEA = dist[fr] - dist[fl]

fl, fr = near(long, 165), near(long, -130)
Length_Ross = dist[fr] - dist[0] + (dist[-1] - dist[fl])

Length_circumpolar = dist[-1] - dist[0]

print("Length A-B = %d km"%Length_AB)
print("Length WAP = %d km"%Length_WAP)
print("Length Weddell = %d km"%Length_Weddell)
print("Length W-EA = %d km"%Length_WEA)
print("Length E-EA = %d km"%Length_EEA)
print("Length Ross = %d km"%Length_Ross)
print("Length circumpolar = %d km"%Length_circumpolar)
