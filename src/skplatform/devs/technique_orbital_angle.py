from typing import Union
import numpy as np
from ..platform import Platform, OpticalGeometry
from ..satellite import SatelliteBase
from ..geodesy import Geodetic


# ------------------------------------------------------------------------------
#           _technique_set_platform_position_from_llh
# ------------------------------------------------------------------------------
def _technique_set_platform_position_from_orbit_plane_angle(platform: Platform,
                                                            ut: np.datetime64,
                                                            position_definition: np.ndarray) -> bool:
    latitude = position_definition[0]
    longitude = position_definition[1]
    height = position_definition[2]
    observer_bearing = position_definition[3]
    if (observer_bearing > 180.0):
        observer_bearing -= 360.0
    observer_height = position_definition[4]
    ok = (latitude > -91.0) and (latitude < 91.0) and \
         (longitude > -361) and (longitude < 721) and \
         (observer_bearing > -361) and (observer_bearing < 721) and \
         (observer_height >= height)
    if (not ok):
        raise ValueError('_technique_set_platform_position_from_observer_looking_at_llh parameters that were out of acceptable range. You might want to check you havethe correct parameter order')

    satellite: SatelliteBase = platform.platform_locator
    ok = isinstance(satellite, SatelliteBase)

    T = satellite.period()
    Tequator = satellite.equator_crossing()
    satellite.eciposition()

    Get orbital anomaly at start time. This is zero
    Get orbital anomaly for one period

    Get eci position at a time
    convert position to orbit angle in range 0 to 360 or higher

    lat0 = zero_lat
    lon0 = np.interp(zero_lat, latitude, longitude)
    geo = Geodetic()
    location = geo.xyz_from_llh([lat0, lon0, 0.0])
    xyz0 = location / np.linalg.norm(location)
    angle_along_orbit = np.ones_like(latitude)
    for idx, (lat, lon) in enumerate(zip(latitude, longitude)):
        location = geo.xyz_from_llh([lat, lon, 0.0])
        xyz = location / np.linalg.norm(location)
        if lat < lat0:
            angle_along_orbit[idx] = -np.arccos(np.dot(xyz, xyz0))
        else:
            angle_along_orbit[idx] = np.arccos(np.dot(xyz, xyz0))

    bearingfunc = ObserverBearingFunction(latitude, longitude, height, observer_height, observer_bearing)               # Initialize the observer bearing calculator function
    minval = 9999.0
    minangle = 9999.0
    for angle in range(-180, 190, 10):                                                                                  # Coarse step through azimuths at the tangent point
        val = math.fabs(bearingfunc(angle))                                                                             # to find the azimuths at the tangent point
        if val < minval:                                                                                                # that are generate the proper bearing required at the observer
            minval = val                                                                                                # the calculation is usually within 10 centimeters.
            minangle = angle

    roots = scipy.optimize.brentq(bearingfunc, minangle - 10, minangle + 10, xtol=1.0E-04, rtol=1.0E-06, full_output=True)
    angle = roots[0]
    converged = roots[1].converged
    if (converged):
        zero = bearingfunc(angle)
        obs = bearingfunc.obslocation
        bearingfunc._geo.from_xyz(obs)
        latitude = bearingfunc._geo.latitude                                                                            # Get the latitude of the observer location
        longitude = bearingfunc._geo.longitude                                                                          # Get the longitude of the observer position
        platform.platform_pointing.set_platform_location(latlonheightandt=(latitude, longitude, observer_height, ut))            # But get the height from the user defined value so it is bang on.
    if not converged:
        logging.warning('technique_set_platform_position_from_observer_looking_at_llh, cannot find observer position that matches the bearing requiremenets')
    return converged

    def latitude(self, orbit_angle):
        lat, lon, angles = self._calipso_position()
        return np.interp(orbit_angle, angles, lat)

    def longitude(self, orbit_angle):
        lat, lon, angles = self._calipso_position()
        return np.interp(orbit_angle, angles, lon)

    def time(self, orbit_angle):
        mjd = self.mjd(orbit_angle)
        return pd.Timedelta(mjd, 'D') + pd.Timestamp('1858-11-17')

    def mjd(self, orbit_angle):
        calipso = xr.open_dataset(self._file, group='CALIPSO')
        lat, lon, angles = self._calipso_position()
        time = calipso.time.values

        mjds = (time - np.datetime64('1858-11-17')) / np.timedelta64(1, 'D')
        mjd = np.interp(orbit_angle, angles, mjds)
        return mjd

    def _calipso_position(self):
        calipso = xr.open_dataset(self._file, group='CALIPSO')
        latitude = np.asarray(calipso.latitude.values, dtype=float)
        longitude = np.asarray(calipso.longitude.values, dtype=float)
        angles = self._angle_along_orbit(latitude, longitude, 0.0)
        return latitude, longitude, angles


def optical_geometry(orbit_angle, altitude=10.0, satellite_altitude=500.0):

        lat = self.latitude(orbit_angle)
        lon = self.longitude(orbit_angle)
        time = self.time(orbit_angle)
        mjd = (time - pd.Timestamp('1858-11-17')) / pd.Timedelta(1, 'D')
        tp = Geodetic()
        location = tp.xyz_from_llh([lat, lon, altitude * 1000])
        normal = self.normal_vector

        look = np.cross(location / np.linalg.norm(location), normal)
        re = 6372
        tp_to_sat = np.sqrt((re + satellite_altitude) ** 2 - (re + altitude) ** 2) * 1000
        obs = tp - look * tp_to_sat
        hor = np.cross(look, obs / np.linalg.norm(obs))
        up = np.cross(hor, look)

        return OpticalGeometry(observer=obs, look_vector=look, local_up=up, mjd=mjd)
