from skorbit.visualization.plotlib import plot_dataset, plot_files, plot_observation, blank_plot

import numpy as np
import xarray as xr
from pathlib import Path

# plot a single file
ds = xr.load_dataset('output/ticfire_2019-08-01_00_01_00_to_2019-08-01_03_00_01.nc')
plot_dataset(ds, 'output/plot_dataset_ticfire_example.png')

# plot only some of the data
plot_dataset(ds, 'output/plot_dataset_ticfire_example2.png', time_step=3, along_step=1, across_step=1)

# Can control the exents so auto-focuses
plot_dataset(ds, 'output/plot_dataset_ticfire_example3.png', time_step=3, along_step=4, across_step=3, extent=None)

# can plot a list of files
plot_files(Path('output').glob('ticfire_*'), 'output/plot_files_ticfire_example.png')

# can plot single observations
time = str(np.datetime_as_string(ds.time[4].values))
plot_observation(ds, 'output/plot_observation_ticfire_example.png', time=time)

# or a list of observations
time = [x for x in range(4,8)]
plot_observation(ds, 'output/plot_observation_ticfire_example3.png', time=time)

time = [str(t) for t in np.datetime_as_string(ds.time[4:8].values)]
plot_observation(ds, 'output/plot_observation_ticfire_example2.png', time=time)

# the same examples but with limb geometry
ds = xr.load_dataset('output/aos-sky_2019-08-01_01_43_24_to_2019-08-01_02_00_27.nc')
plot_dataset(ds, 'output/plot_dataset_aos-sky_example.png')

plot_dataset(ds, 'output/plot_dataset_aos-sky_example2.png', time_step=60, along_step=4, across_step=3)

ds = xr.load_dataset('output/aos-sky_2019-08-01_02_00_28_to_2019-08-01_02_17_31.nc')
plot_dataset(ds, 'output/plot_dataset_aos-sky_example3.png', time_step=60, along_step=4, across_step=3, extent=None)


plot_files(Path('output').glob('aos-sky_*'), 'output/plot_files_aos-sky_example.png')

time = '2019-08-01 00:51:00'
plot_observation(ds, 'output/plot_observation_aos-sky_example.png', time=time)

# can pass in an existing axes object to combine plots for things like co-location 
aossky_ds = xr.load_dataset('output/single_file_a.nc')
aossky_time = '2019-08-01 00:06:15'
ticfire_ds = xr.load_dataset('output/ticfire_2019-08-01_00_01_00_to_2019-08-01_03_00_01.nc')
ticfire_time = '2019-08-01 00:01:00'

fig, ax = blank_plot()
plot_observation(ticfire_ds, '', time=ticfire_time, axes=ax)
# if passed an axes object, the image will not be saved
plot_observation(aossky_ds, 'output/this_plot_will_not_be_made.png', time=aossky_time, axes=ax)
fig.savefig('output/combined_observation_example.png')


aossky2_ds = xr.load_dataset('output/imaginary_sky.nc')
plot_dataset(aossky2_ds, 'output/plot_aossky2_ds_example.png', extent=None)
