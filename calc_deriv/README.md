## Codes for calculating derived variables from POP/CICE outputs

This directory contains codes that were used to process the raw (monthly- and daily- averaged) POP/CICE outputs. Most of them produce one or more of the \*.npz and netCDF files located in the `../data_reproduce_figs` directory.

Brief description of each \*.py file:

* `201a-avg_vortbdgt.py`: Calculates seasonal averages and the 2005-2009 average of the depth-dependent vorticity balance terms.
* `201b-int_avg_vortbdgt.py`: Vertically-integrates the time-averaged, depth-dependent vorticity balance terms (*i.e.*, the outputs of the above script).
* `201e-calc_vortbdgt_daily.py`: Vertically-integrates the daily, depth-dependent vorticity balance terms.
* `206a-calc_circulation_circumpolar.py`: Calculates a time series of vertically-integrated vorticity balance terms (*i.e.*, the outputs from the above script) area-averaged within the circumpolar control volume formed by the model 800 m and 2500 m isobaths.
* `206b-calc_circulation_segments.py`: Calculates the same time series as the above script, but area-averaged within segments of the circumpolar control volume (Amundsen-Bellingshausen, Western Antarctic Peninsula, Weddell, Western East Antarctica, Eastern East Antarctica and Ross, see Figure 1 in the manuscript).
* `210a-calc_wind-ice-ocean_stress.py`: Calculates the surface ocean stress vector (stress due to the wind plus stress due to the relative sea ice motion) from POP and CICE outputs.
* `300-calc_zIT.py`: Calculates the vertical index of the deepest cell in the POP model grid that has no sidewalls.

Please contact [André Palóczy](mailto:apaloczy@ucsd.edu) with any questions on the code.
