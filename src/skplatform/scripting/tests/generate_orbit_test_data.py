from typing import Union, Tuple

import sys
import numpy as np
from numpy.typing import NDArray
from pathlib import Path
from skplatform import Platform
from skplatform.satellite import SatelliteSGP4
from datetime import timezone

from skplatform.scripting.rotations import calculate_2d_rotation_matrices
from skplatform.scripting.tests.utils import log10_histogram

"""
This script is intended to generate the data needed for test_orbits.py. The values 
saved in the text files under ./tests should be used to update the arrays stored in 
orbit_test_data.py. The additional package skyfield must be installed, and the 
existence of the ./tests folder is assumed.
"""


TEST_DIR = Path(__file__).parent


def ferrari(x, y, z):
    """
    Converts x,y,z geographic coordinates to geodetic.
    Simply here as a backup test of the xyz->geodetic conversion

    Only returns the latitude.

    From https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
    """
    a = 6378137.0
    b = 6356752.3142
    e2 = (a*a - b*b) / (a*a)
    e2p = (a*a - b*b) / (b*b)
    p = np.sqrt(x*x + y*y)
    F = 54 * b*b*z*z
    G = p*p + (1-e2)*z*z - e2 * (a*a - b*b)
    c = e2*e2*F*p*p / (G**3)
    s = np.cbrt(1.0 + c + np.sqrt(c*c + 2*c))
    k = s + 1.0 + 1/s
    P = F / (3 * k*k*k*G*G)
    Q = np.sqrt(1 + 2*e2*e2*P)
    r0 = -(P*e2*p) / (1.0+Q) + np.sqrt(0.5*a*a*(1+1/Q) - (P*(1-e2)*z*z)/(Q*(1+Q)) - 0.5*P*p*p)
    U = np.sqrt((p - e2*r0)**2 + z*z)
    V = np.sqrt((p - e2*r0)**2 + (1-e2)*z*z)
    z0 = b*b*z / (a*V)
    h = U * (1.0 - b*b/(a*V))  # altitude
    l = np.arctan2(y, x)  # longitude
    return np.arctan((z + e2p*z0) / p)  # latitude


def setup_sat_from_TLE(tle_file: Path) -> SatelliteSGP4:
    """_summary_
    Returns a Skplatform SatelliteSGP4 instance. The satellite will have polar motion 
    correction enabled.

    Parameters
    ----------
    tle_file : Path
        _description_ The path and filename of the TLE to initialize the satellite with.

    Returns
    ------- 
    _type_ : SatelliteSGP4
        _description_ A satellite instance with polar motion corrections enabled.
    """
    with open(tle_file) as file:
        lines = file.readlines()
    sat = SatelliteSGP4(twolines=lines)
    sat.use_polar_motion = True
    return sat


def create_skyfield_sat():
    """
    Returns a skyfield satellite, and the geodetic package and timescale for using it. 
    The geodetic package and timescale are needed to convert reference frames. 
    """
    from skyfield.api import EarthSatellite, load, wgs84
    from skyfield.data import iers

    url = load.build_url('finals2000A.all')
    with load.open(url) as f:
        finals_data = iers.parse_x_y_dut1_from_finals_all(f)

    ts = load.timescale()
    iers.install_polar_motion_table(ts, finals_data)
    
    # To use the WGS84 gravity model it says to use SGP4 Satrec
    # https://rhodesmill.org/skyfield/earth-satellites.html
    from sgp4.api import Satrec, WGS84

    with open('../../../AOS_Sky_Descending_v2025-01-29_Epoch2019.tle') as file:
        lines = file.readlines()
    
    # sat = EarthSatellite(lines[0], lines[1], 'Test Dummy', ts)
    sgp4_sat = Satrec.twoline2rv(lines[0], lines[1], WGS84)
    sat = EarthSatellite.from_satrec(sgp4_sat, ts)
    
    return sat, wgs84, ts  # need ts for time conversions and wgs84 for geometry


def absolute_difference(a, b):
    """_summary_
    Calculates the largest absolute difference in magnitude of two arrays. 
    The formula is max(|a-b|).

    Parameters
    ----------
    a : _type_ : NDArray
        _description_ A Numpy array, or an object Numpy can coerce into an array.
    b : _type_ : NDArray
        _description_ A Numpy array, or an object Numpy can coerce into an array.

    Returns
    -------
    _type_ : float
        _description_ The largest difference of the two inputs, by magnitude. 
    """
    return np.nanmax(np.abs( (a-b)))


def percent_difference(a, b):
    """_summary_
    Calculates the largest percent difference in magnitude of two arrays. 
    The formula is max(|a-b|/b)*100.0.

    Parameters
    ----------
    a : _type_ : NDArray
        _description_ A Numpy array, or an object Numpy can coerce into an array.
    b : _type_ : NDArray
        _description_ A Numpy array, or an object Numpy can coerce into an array.

    Returns
    -------
    _type_ : float
        _description_ The largest percent difference of the two inputs, by magnitude.
    """
    return np.nanmax(np.abs((a-b) / b)) * 100.0


def print_test_header(name: str, desc: str):
    """
    Generates a 4 line text header composed of a line of = characters, the name, the 
    description and a closing line of = characters. Used to provide visual separation
    of values printed to the console.
    """
    print('')
    print('='*80)
    print(name)
    print(desc)
    print('='*80)


def write_array_to_txt(filename: Tuple[str, Path], arr: NDArray):
    """_summary_
    Writes the passed-in array into a text file, as a numpy array definition, ready
    to be pasted into code. This function is meant to make updating the orbit 
    unit tests easier. The precision used is 32, to minimize floating point error
    while saving and loading the values.

    Parameters
    ----------
    filename : Tuple[str, Path]
        _description_ The filename and path to save the array with.
    arr : NDArray
        _description_ The array to save. It assumed to be a float array.
    """
    with open(filename, 'w') as file: 
        # according to Numpy documentation setting precision to None
        # means the precision required to represent the number will
        # be used but this didn't seem to work, there were still 
        # floating point errors.
        with np.printoptions(precision=32, threshold=sys.maxsize):
            file.write('np.')
            file.write(np.array_repr(arr.flatten()))


def check_eci(scn: Platform, tle_file: Union[Path, str]):
    """
    Compares the ECI position of the satellite calculated using skplatform, Astropy and Skyfield. 
    The computed values are saved in text files to allow updating orbit_test_data.py. 
    """
    from skplatform import Platform
    from sgp4.api import Satrec, WGS84

    from astropy.time import Time
    from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
    from astropy import units as u

    from skyfield.api import Distance, Velocity
    from skyfield.positionlib import Geocentric
    from skyfield.sgp4lib import TEME as skyTEME
    
    sky_sat, _, ts = create_skyfield_sat()
   
    # get the skplatform velocity and positions
    sat = setup_sat_from_TLE(tle_file)
    sk_plat = Platform(platform_locator=sat)

    tle = open(tle_file).readlines()
    astro_sat = Satrec.twoline2rv(tle[0], tle[1], WGS84)

    # store the eci coordinates
    sat_eci = np.zeros((scn.num_measurements, 3))
    astro_eci = np.zeros((scn.num_measurements, 3))
    sky_eci = np.zeros((scn.num_measurements, 3))
    sky_sat_eci = np.zeros((scn.num_measurements, 3))
    astro_sat_eci = np.zeros((scn.num_measurements, 3))

    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]
        # SkPlatform doesn't store the eci values 
        sk_plat.platform_locator.update_eci_position(time)

        sat_eci_pos = sk_plat.platform_locator.eciposition()
        sat_eci_vel = sk_plat.platform_locator.ecivelocity()
        teme_eci_pos = CartesianRepresentation(sat_eci_pos * u.m)
        teme_eci_vel = CartesianDifferential(sat_eci_vel * u.m/u.s)
        
        astro_time = Time(time)
        teme = TEME(teme_eci_pos.with_differentials(teme_eci_vel), obstime=astro_time)

        # get the astropy sgp4 values
        err, astro_pos, astro_v = astro_sat.sgp4(astro_time.jd1, astro_time.jd2)  # km and km/s

        # sky_geocentric is in GCRS internally, so everything in/out needs conversion
        sky_time = ts.utc(time.replace(tzinfo=timezone.utc))
        sky_pos = Distance(m=sat_eci_pos)
        sky_vel = Velocity.m_per_s(sat_eci_vel)
        sky_geocentric = Geocentric.from_time_and_frame_vectors(sky_time, skyTEME, sky_pos, sky_vel)

        sat_eci[t, :] = sat_eci_pos
        sky_eci[t, :] = sky_geocentric.frame_xyz(skyTEME).m  # internally it is in GCRS, so we need to convert back
        sky_sat_eci[t, :] = sky_sat._position_and_velocity_TEME_km(sky_time)[0]*1000
        astro_eci[t, :] = teme.cartesian.xyz.to_value()
        astro_sat_eci[t, :] = np.array(astro_pos)*1000.0

    print_test_header('check_eci()', 'Checking calculation of ECI coordinates')
    print(f'Comparison of ECI (TEME) Coordinates\n'
          f'\tskplatform vs Skyfield:         {absolute_difference(sat_eci.flatten(), sky_eci.flatten()):.5e} m'
          f'\t{percent_difference(sat_eci.flatten(), sky_eci.flatten()):.5e}%\n'
          f'\tskplatform vs sk->Astropy:      {absolute_difference(sat_eci.flatten(), astro_eci.flatten()):.5e} m'
          f'\t{percent_difference(sat_eci.flatten(), astro_eci.flatten()):.5e}%\n'
          f'\tskplatform vs Skyfield Sat:     {absolute_difference(sat_eci.flatten(), sky_sat_eci.flatten()):.5e} m'
          f'\t{percent_difference(sat_eci.flatten(), sky_sat_eci.flatten()):.5e}%\n'
          f'\tskplatform vs SGP4 Sat:         {absolute_difference(sat_eci.flatten(), astro_sat_eci.flatten()):.5e} m'
          f'\t{percent_difference(sat_eci.flatten(), astro_sat_eci.flatten()):.5e}%\n'
          )
    
    write_array_to_txt(TEST_DIR.joinpath('data/sat_eci.txt'), sat_eci)
    write_array_to_txt(TEST_DIR.joinpath('data/sky_eci.txt'), sky_eci)
    write_array_to_txt(TEST_DIR.joinpath('data/astro_eci.txt'), astro_eci)
    write_array_to_txt(TEST_DIR.joinpath('data/sky_sat_eci.txt'), sky_sat_eci)
    write_array_to_txt(TEST_DIR.joinpath('data/astro_sat_eci.txt'), astro_sat_eci)
    

def check_eci_to_itrs(scn: Platform, tle_file: Union[Path, str]):
    """
    Compares the ECEF (ITRS) position of the satellite calculated using skplatform, Astropy and Skyfield. 
    The computed values are saved in text files to allow updating orbit_test_data.py.
    """
    from skplatform import Platform

    from astropy.time import Time
    from astropy import units as u
    from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation, ITRS
    
    from skyfield.api import Distance, Velocity
    from skyfield.positionlib import Geocentric
    from skyfield.framelib import itrs
    from skyfield.sgp4lib import TEME as skyTEME
    
    # only need the timescale object
    _, _, ts = create_skyfield_sat()
   
    # get the skplatform velocity and positions
    sat = setup_sat_from_TLE(tle_file)
    sk_plat = Platform(platform_locator=sat)

    sat_ecef = np.zeros((scn.num_measurements, 3))
    astro_ecef = np.zeros((scn.num_measurements, 3))
    sky_ecef = np.zeros((scn.num_measurements, 3))

    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]
        # SkPlatform doesn't store the eci values 
        sk_plat.platform_locator.update_eci_position(time)

        sat_eci_pos = sk_plat.platform_locator.eciposition()
        sat_eci_vel = sk_plat.platform_locator.ecivelocity()
        teme_eci_pos = CartesianRepresentation(sat_eci_pos * u.m)
        teme_eci_vel = CartesianDifferential(sat_eci_vel * u.m/u.s)
        
        astro_time = Time(time)

        teme = TEME(teme_eci_pos.with_differentials(teme_eci_vel), obstime=astro_time)
        itrs_geo = teme.transform_to(ITRS(obstime=astro_time))

        # sky_geocentric is in GCRS internally, so everything in/out needs conversion
        sky_time = ts.utc(time.replace(tzinfo=timezone.utc))
        sky_pos = Distance(m=sat_eci_pos)
        sky_vel = Velocity.m_per_s(sat_eci_vel)
        sky_geocentric = Geocentric.from_time_and_frame_vectors(sky_time, skyTEME, sky_pos, sky_vel)

        sat_ecef[t, :] = sk_plat.platform_locator.position
        sky_ecef[t, :] = sky_geocentric.frame_xyz(itrs).m
        astro_ecef[t, :] = itrs_geo.cartesian.xyz.to_value()

    print_test_header('check_eci_to_itrs()', 'Checking ECI to ITRS conversion')
    
    print(f'Comparison of ITRS Coordinates\n'
          f'\tskplatform vs Skyfield: {absolute_difference(sat_ecef, sky_ecef):.5e} m\t{percent_difference(sat_ecef, sky_ecef):.5e}%\n'
          f'\tskplatform vs Astropy:  {absolute_difference(sat_ecef, astro_ecef):.5e} m\t{percent_difference(sat_ecef, astro_ecef):.5e}%\n'
          )
    
    write_array_to_txt(TEST_DIR.joinpath('data/sat_ecef.txt'), sat_ecef)
    write_array_to_txt(TEST_DIR.joinpath('data/sky_ecef.txt'), sky_ecef)
    write_array_to_txt(TEST_DIR.joinpath('data/astro_ecef.txt'), astro_ecef)


def check_itrs_to_lla(scn: Platform, tle_file: Union[Path, str]):
    """
    Compares the latitude, longitude and altitude of the satellite calculated using skplatform, Astropy and Skyfield. 
    The input is the ECEF (ITRS) position calculated by skplatform.
    The computed values are saved in text files to allow updating orbit_test_data.py.
    """
    from skplatform import Platform
    
    from astropy.coordinates import CartesianRepresentation
    from astropy import units as u
    from astropy.coordinates import ITRS

    from skyfield.api import Distance, Velocity
    from skyfield.positionlib import Geocentric
    from skyfield.framelib import itrs

    _, wgs84, ts = create_skyfield_sat()
   
    # get the skplatform velocity and positions
    sat = setup_sat_from_TLE(tle_file)
    sk_plat = Platform(platform_locator=sat)

    sk_lla = np.zeros((scn.num_measurements, 3))
    astro_lla = np.zeros((scn.num_measurements, 3))
    sky_lla = np.zeros((scn.num_measurements, 3))
    
    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]
        # SkPlatform doesn't store the eci values 
        sk_plat.platform_locator.update_eci_position(time)

        astro_itrs = ITRS(CartesianRepresentation(sk_plat.platform_locator.position * u.m))
        astro_geo = astro_itrs.earth_location.geodetic
        astro_lla[t, :] = np.array([astro_geo.lat.value, astro_geo.lon.value,  astro_geo.height.value])
        
        sky_time = ts.utc(time.replace(tzinfo=timezone.utc))
        sky_pos = Distance(m=sk_plat.platform_locator.position)
        sky_vel = Velocity.m_per_s(sk_plat.platform_locator.velocity)
        sky_geocentric = Geocentric.from_time_and_frame_vectors(sky_time, itrs, sky_pos, sky_vel)
        sky_geodetic = wgs84.geographic_position_of(sky_geocentric)
        sky_lla[t, :] = np.array([sky_geodetic.latitude.degrees, sky_geodetic.longitude.degrees, sky_geodetic.elevation.m])
    
    scn.propagate()
    sk_lla[:, 0] = scn.observer_latitude
    sk_lla[:, 1] = scn.observer_longitude
    sk_lla[:, 2] = scn.observer_altitude

    print_test_header('check_itrs_to_lla()', 'Checking ECEF to Geodetic conversion')
    print(f'skplatform vs Astropy:\n'
          f'\tLat: {absolute_difference(sk_lla[:,0].flatten(), astro_lla[:,0].flatten()):.5e} m\t{percent_difference(sk_lla[:,0].flatten(), astro_lla[:,0].flatten()):.5e}%\n'
          f'\tLon: {absolute_difference(sk_lla[:,1].flatten(), astro_lla[:,1].flatten()):.5e} m\t{percent_difference(sk_lla[:,1].flatten(), astro_lla[:,1].flatten()):.5e}%\n'
          f'\tAlt: {absolute_difference(sk_lla[:,2].flatten(), astro_lla[:,2].flatten()):.5e} m\t{percent_difference(sk_lla[:,2].flatten(), astro_lla[:,2].flatten()):.5e}%\n'
          f'skplatform vs Skyfield:\n'
          f'\tLat: {absolute_difference(sk_lla[:,0].flatten(), sky_lla[:,0].flatten()):.5e} m\t{percent_difference(sk_lla[:,0].flatten(), sky_lla[:,0].flatten()):.5e}%\n'
          f'\tLon: {absolute_difference(sk_lla[:,1].flatten(), sky_lla[:,1].flatten()):.5e} m\t{percent_difference(sk_lla[:,1].flatten(), sky_lla[:,1].flatten()):.5e}%\n'
          f'\tAlt: {absolute_difference(sk_lla[:,2].flatten(), sky_lla[:,2].flatten()):.5e} m\t{percent_difference(sk_lla[:,2].flatten(), sky_lla[:,2].flatten()):.5e}%\n'
          )

    write_array_to_txt(TEST_DIR.joinpath('data/sk_lla.txt'), sk_lla)
    write_array_to_txt(TEST_DIR.joinpath('data/astro_lla.txt'), astro_lla)
    write_array_to_txt(TEST_DIR.joinpath('data/sky_lla.txt'), sky_lla)


def check_orbit_to_lla(scn: Platform, tle_file: Union[Path, str]):
    """
    Compares the latitude, longitude and altitude of the satellite calculated using skplatform, Astropy and Skyfield. 
    The values are calculated from the observation times using each library's propagator interface. The results are
    the other libraries are independent of skplatform.

    Compares skplatform geodetic to Astropy Geodetic and skplatform Geodetic vs Sk ECI -> TEME -> ITRS -> Geodetic
    
    The computed values are saved in text files to allow updating orbit_test_data.py.
    """
    from datetime import timezone
    from astropy.time import Time
    from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
    from astropy import units as u
    from astropy.coordinates import ITRS
    from sgp4.api import Satrec, WGS84

    sky_sat, wgs84, ts = create_skyfield_sat()

    tle = open(tle_file).readlines()
    astro_sat = Satrec.twoline2rv(tle[0], tle[1], WGS84)

    scn.propagate()

    astro_geo = np.zeros((scn.num_measurements, 3))  # lat, lon, alt
    sky_geo = np.zeros((scn.num_measurements, 3))

    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]

        astro_time = Time(time.replace(tzinfo=timezone.utc))
        err, astro_p, astro_v = astro_sat.sgp4(astro_time.jd1, astro_time.jd2)  # ECI, km and km/s

        astro_pos = CartesianRepresentation(astro_p * u.km)
        astro_vel = CartesianDifferential(astro_v * u.km/u.s)
        astro_teme = TEME(astro_pos.with_differentials(astro_vel), obstime=astro_time)
        astro_loc = astro_teme.transform_to(ITRS(obstime=astro_time)).earth_location
        
        astro_geo[t, 0] = astro_loc.geodetic.lat.to_value()
        astro_geo[t, 1] = astro_loc.geodetic.lon.to_value()
        astro_geo[t, 2] = astro_loc.geodetic.height.to(u.m).to_value()

        sky_time = ts.utc(time.replace(tzinfo=timezone.utc))
        sky_geocentric = sky_sat.at(sky_time)
        sky_geodetic = wgs84.geographic_position_of(sky_geocentric)
        sky_geo[t, :] = np.array([sky_geodetic.latitude.degrees, sky_geodetic.longitude.degrees, sky_geodetic.elevation.m]) 

    print_test_header('check_orbit_to_lla()',
                      'Compare Geodetic coordinates directly from propagators')

    print(f'Skplatform vs SGP4->Astropy\n'
          f'\tLat: {absolute_difference(scn.observer_latitude, astro_geo[:,0].flatten()):.5e} deg\t{percent_difference(scn.observer_latitude, astro_geo[:,0].flatten()):.5e}%\n'
          f'\tLon: {absolute_difference(scn.observer_longitude, astro_geo[:,1].flatten()):.5e} deg\t{percent_difference(scn.observer_longitude, astro_geo[:,1].flatten()):.5e}%\n'
          f'\tAlt: {absolute_difference(scn.observer_altitude, astro_geo[:,2].flatten()):.5e} m\t{percent_difference(scn.observer_altitude, astro_geo[:,2].flatten()):.5e}%\n'
          f'Skplatform vs Skyfield\n'
          f'\tLat: {absolute_difference(scn.observer_latitude, sky_geo[:,0].flatten()):.5e} deg\t{percent_difference(scn.observer_latitude, sky_geo[:,0].flatten()):.5e}%\n'
          f'\tLon: {absolute_difference(scn.observer_longitude, sky_geo[:,1].flatten()):.5e} deg\t{percent_difference(scn.observer_longitude, sky_geo[:,1].flatten()):.5e}%\n'
          f'\tAlt: {absolute_difference(scn.observer_altitude, sky_geo[:,2].flatten()):.5e} m\t{percent_difference(scn.observer_altitude, sky_geo[:,2].flatten()):.5e}%\n'
          )

    write_array_to_txt(TEST_DIR.joinpath('data/sk_lat.txt'), scn.observer_latitude)
    write_array_to_txt(TEST_DIR.joinpath('data/sk_lon.txt'), scn.observer_longitude)
    write_array_to_txt(TEST_DIR.joinpath('data/sk_alt.txt'), scn.observer_altitude)
    write_array_to_txt(TEST_DIR.joinpath('data/astro_lat.txt'), astro_geo[:, 0])
    write_array_to_txt(TEST_DIR.joinpath('data/astro_lon.txt'), astro_geo[:, 1])
    write_array_to_txt(TEST_DIR.joinpath('data/astro_alt.txt'), astro_geo[:, 2])
    write_array_to_txt(TEST_DIR.joinpath('data/sky_lat.txt'), sky_geo[:, 0])
    write_array_to_txt(TEST_DIR.joinpath('data/sky_lon.txt'), sky_geo[:, 1])
    write_array_to_txt(TEST_DIR.joinpath('data/sky_alt.txt'), sky_geo[:, 2])


def check_orbit_is_planar(scn: Platform, tle_file: Union[Path, str]):
    """
    Test how planar the orbit is. An ideal orbit would reside entirely on a single plane when plotted in the Earth 
    Centered Inertial (ECI) frame. The SGP4 propagator includes perturbations that prevent the orbit from being 
    perfectly planar, but a significant deviance is not expected. Reference values were not found, but values of 
    0.3 degrees deviation were calculated. A perfectly planar orbit or an orbit with significant deviance indicate 
    a problem. 

    Planarity is checked by comparing the angle between the cross-product of position and velocity of each observation 
    to that of every other observation.

    The computed values are saved in text files to allow updating orbit_test_data.py.
    """
    from skplatform import Platform

    sat = setup_sat_from_TLE(tle_file)
    teme_plat = Platform(platform_locator=sat)
    cross_vectors = np.zeros((scn.num_measurements, 3))
    for time in range(scn.num_measurements):
        # SkPlatform doesn't store these in an array
        teme_plat.platform_locator.update_eci_position(scn.measurement_times[time])
        teme_pos = teme_plat.platform_locator.eciposition()
        teme_pos /= np.linalg.norm(teme_pos)
        teme_vel = teme_plat.platform_locator.ecivelocity()
        teme_vel /= np.linalg.norm(teme_vel)
        axis = np.cross(teme_pos, teme_vel)
        cross_vectors[time, :] = axis / np.linalg.norm(axis)

    compare = []
    for i in range(scn.num_measurements):
        for j in range(i+1, scn.num_measurements):
            compare.append(1.0 - np.dot(cross_vectors[i], cross_vectors[j]))

    print_test_header('check_orbit_is_planar()', 
                      'Check for planar orbit in skplatform ECI frame (compare position cross '
                      'velocity between all measurements)')
    print(f'Difference {np.max(np.abs(compare))*100:.5e}'
          f'\t{np.rad2deg(np.arccos(1.0 - np.max(compare))):.5e} deg'
          f'\tMean: {np.mean(compare)*100.0:.5e}'
          f'\tMedian: {np.median(compare)*100.0:.5e}')

    print('Histogram:')
    normal_vector_histogram, histogram_header = log10_histogram(np.array(compare), header=True, 
                                                                width=5, smallest_power=-16, num_bins=13)
    print(histogram_header, 'Log10')
    print(normal_vector_histogram, 'Position cross Velocity Difference')
    print('')

    write_array_to_txt(TEST_DIR.joinpath('data/planar.txt'), np.array(compare))


def setup_limb_propagator(tle_file: Union[Path, str]):
    """_summary_
    Creates an instance of the LimbPropagator for testing.

    Parameters
    ----------
    tle_file : Union[Path, str]
        _description_ The filename and path of the TLE to use ot initialize the satellite.

    Returns
    -------
    _type_ : LimbPropagator
        _description_ A fully initialized LimbPropagator instance ready to be used in testing.
    """
    from datetime import timedelta
    from skplatform import Platform
    from skplatform.scripting.propagator import LimbPropagator

    sat = setup_sat_from_TLE(tle_file)
    platform = Platform(platform_locator=sat)
    # approximately 5 orbits, with measurements every 100s
    times = np.array([timedelta(seconds=dt) for dt in range(0, 90*60*5, 100)])
    rot_mat = calculate_2d_rotation_matrices(0, 4.9656*2, 1, 39)
    boresight = 180.0
    altitude = 20_000.0
    limb_propagator = LimbPropagator(platform, times, rot_mat, boresight, altitude)
    return limb_propagator


def setup_nadir_propagator(tle_file: Union[Path, str]):
    """_summary_
    Creates an instance of the NadirPropagator for testing.

    Parameters
    ----------
    tle_file : Union[Path, str]
        _description_ The filename and path of the TLE to use ot initialize the satellite.

    Returns
    -------
    _type_ : NadirPropagator
        _description_ A fully initialized NadirPropagator instance ready to be used in testing.
    """
    from datetime import timedelta
    from skplatform import Platform
    from skplatform.scripting.propagator import NadirPropagator

    sat = setup_sat_from_TLE(tle_file)
    platform = Platform(platform_locator=sat)
    # approximately 5 orbits, with measurements every 100s
    times = np.array([timedelta(seconds=dt) for dt in range(0, 90*60*5, 100)])
    rot_mat = calculate_2d_rotation_matrices(0, 4.9656 * 2, 1, 39)
    nadir_propagator = NadirPropagator(platform, times, rot_mat)

    return nadir_propagator


def check_limb_orbit_conversions():
    """
    Runs the frame conversion comparisons using the LimbPropagator class.
    """
    tle_file = '../../../AOS_Sky_Descending_v2025-01-29_Epoch2019.tle'
    limb_propagator = setup_limb_propagator(tle_file)
    check_eci(limb_propagator, tle_file)
    limb_propagator = setup_limb_propagator(tle_file)
    check_eci_to_itrs(limb_propagator, tle_file)
    limb_propagator = setup_limb_propagator(tle_file)
    check_itrs_to_lla(limb_propagator, tle_file)
    limb_propagator = setup_limb_propagator(tle_file)
    check_orbit_to_lla(limb_propagator, tle_file)


def check_nadir_orbit_conversions():
    """
    Runs the frame conversion comparisons using the NadirPropagator class.
    """
    tle_file = '../../../AOS_Sky_Descending_v2025-01-29_Epoch2019.tle'
    nadir_propagator = setup_nadir_propagator(tle_file)
    check_eci(nadir_propagator, tle_file)
    nadir_propagator = setup_nadir_propagator(tle_file)
    check_eci_to_itrs(nadir_propagator, tle_file)
    nadir_propagator = setup_nadir_propagator(tle_file)
    check_itrs_to_lla(nadir_propagator, tle_file)
    nadir_propagator = setup_nadir_propagator(tle_file)
    check_orbit_to_lla(nadir_propagator, tle_file)


def check_limb_orbit_is_planar():
    """
    Check that the LimbPropagator orbit is near-planar
    """
    tle_file = '../../../AOS_Sky_Descending_v2025-01-29_Epoch2019.tle'
    limb_propagator = setup_limb_propagator(tle_file)
    check_orbit_is_planar(limb_propagator, tle_file)


def check_nadir_orbit_is_planar():
    """
    Check that the NadirPropagator orbit is near-planar
    """
    tle_file = '../../../AOS_Sky_Descending_v2025-01-29_Epoch2019.tle'
    nadir_propagator = setup_nadir_propagator(tle_file)
    check_orbit_is_planar(nadir_propagator, tle_file)


def check_limb_orbit():
    """
    Runs all the checks on the LimbPropagator.
    """
    check_limb_orbit_conversions()
    check_limb_orbit_is_planar()


def check_nadir_orbit():
    """
    Runs all the checks on the NadirPropagator.
    """
    check_nadir_orbit_conversions()
    check_nadir_orbit_is_planar()


if __name__ == "__main__":
    print('Comparing the Nadir and Limb propagator classes. There should be no difference between the two.')
    check_nadir_orbit()
    check_limb_orbit()
