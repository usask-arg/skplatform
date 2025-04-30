from skplatform.scripting.interface import simulate_orbit, open_tle
from skplatform.scripting.conversion import convert_to_single_file, convert_to_single_orbit_files
from skplatform.scripting.instruments import LimbInstrument

import numpy as np
from pathlib import Path
from datetime import datetime

EXAMPLE_DIR = Path(__file__).parent

# output = Path('this_is_a_test.nc')
output = EXAMPLE_DIR.joinpath('output')  # save to the output directory, this is relative to where the script is being ran from
# use a short time period
start = datetime(year=2019, month=8, day=1, hour=0, minute=1, second=0)
end = datetime(year=2019, month=8, day=1, hour=0, minute=5, second=0)


# this line is required for multi-process, otherwise it crashes
if __name__ == '__main__':
    # demonstrate loading a tle file and the creation of a custom instrument
    tle = EXAMPLE_DIR.joinpath('AOS_Storm_Inclined_v2025-01-29_Epoch2019.tle')
    tle_lines = open_tle(tle)
    # The creation of instruments with custom properties is supported by the LimbInstrument and
    # NadirInstrument classes.

    # This defines a limb instrument with a 0.5 degree FOV along the satellite track (ie. vertical FOV) and
    # a 7 degree FOV across the track (horizontal). The simulator will divide the along track FOV into 5
    # angles. Both endpoints are included. The range extends from +FOV/2 to -FOV/2, so the angles for this
    # instrument will be [0.25, 0.125, 0.0, -0.125, -0.25]. If an even number of pixels is requested then
    # there will be no measurement(s) corresponding to the central boresight.

    # Since the frequency argument is an array the simulator will use a non-uniform interval between
    # observations. The first observation will be calculated, 2.0 seconds later a 2nd observation will be
    # calculated, 5.0 seconds after the second observation the third will take place. Then there will be a
    # 9.0 second interval before the cycle begins again.

    # The look_angle_deg argument is the angle between the satellite's velocity, and the central boresight.
    # The target_altitude_m is set to 20 km. So the central boresight is backward looking at 20 km.

    # Since output is set to None and a single thread is being used, this call will only return the
    # dataset, it will not be saved.
    inst = LimbInstrument('inst1', tle_lines, along_fov=0.50, across_fov=7.0, along_pixels=5, across_pixels=30,
                               frequency=np.array([2.0, 5.0, 9.0]), look_angle_deg=180.0, target_altitude_m=20_000.0)
    ds = simulate_orbit(instrument=inst, output=None, tle=tle, start=start, end=end, threads=1)

    # single threaded ticfire example, with conversion of output into orbit files.
    ds = simulate_orbit(instrument='ticfire', output=None, tle=tle, start=start, end=end, threads=1)
    convert_to_single_file(output, output.joinpath('single_file_t.nc'), 'ticfire_*')
    convert_to_single_orbit_files(output, output, 'ticfire_*', output_prefix='orbit_t')

    # multi-threaded AOS-Sky example
    # simulate_orbit() will output files named 'aos-sky_{start}_to_{end}.nc
    ds = simulate_orbit(instrument='aos-sky', output=output, tle=tle, start=start, end=end, threads=8)
    # in multi-threaded mode sd will be None
    print(f'The returned dataset is: {ds}')

    # Make another 2D detector, but only make observations every minute
    tle_lines = open_tle(tle)
    inst2 = LimbInstrument('inst2', tle_lines,
                           along_fov=0.50, across_fov=7.0,
                           along_pixels=5, across_pixels=30,
                           frequency=60.0,
                           look_angle_deg=180.0,
                           target_altitude_m=20_000.0
                           )

    # use the default start time, and run for 3276.8 minutes
    # this will generate enough for 3 full threads and a partial
    # four files should be generated
    simulate_orbit(instrument=inst2, output=output,
                   tle=tle, length=3.2*1024*60.0, threads=4)
    # compact all four files into one
    convert_to_single_file(output, output.joinpath('i2_single_file.nc'), 'inst2_*')
    # separate the files into individual orbits
    convert_to_single_orbit_files(output, output, 'inst2_*', output_prefix='i2_orbit')
