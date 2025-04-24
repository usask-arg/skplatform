from typing import Tuple, List, Optional, Any
from numpy.typing import NDArray

from abc import ABC, abstractmethod
import numpy as np
import xarray as xr
from xarray import Dataset
from datetime import datetime, timedelta

from skplatform.platform import Platform
from skplatform.geodesy import Geodetic

import logging


# import traceback
# import warnings
# import sys

# def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

#     log = file if hasattr(file,'write') else sys.stderr
#     traceback.print_stack(file=log)
#     log.write(warnings.formatwarning(message, category, filename, lineno, line))

# warnings.showwarning = warn_with_traceback

''' Propagator shouldn't be responsible for general rotation matrices, just their application'''

'''
if we are ever required to develop solutions for pointing, consider adding a matrix in the ICF frame
to do the adjustment. Make this an array or a callable so that can incorporate things like jitter
'''

class Propagator(ABC):
    def __init__(self, platform: Platform, time_offsets: NDArray[timedelta], fov_rotation_matrix: NDArray):
        """
        The Propagator subclasses manages propagation of the orbit,
        calculation of the lines of sight, storage the data
        generated and manages the detector geometry.

        The Propagator class is a base class not meant to be used directly.

        Subclasses of the Propagator class are to be created by the
        Instrument classes.
        """
        self.platform: Platform = platform
        self.base_time: datetime = platform.platform_locator.time
        self.measurement_times: NDArray[datetime] = self.base_time + time_offsets
        self.num_measurements: int = len(self.measurement_times)
        # TODO: need to check that R is [along, across, 3, 3]
        self.R: NDArray[float] = fov_rotation_matrix
        self.num_along: int = self.R.shape[0]
        self.num_across: int = self.R.shape[1]
        self.ds : Optional[Dataset] = None

    @abstractmethod
    def propagate():
        """
        Subclasses must implement this function. `propgate()` Should
        propagate the orbit, calculate the location of measurements
        and store the results back into class variables. Should
        call self.to_dataset() to set self.ds to the results of the
        propagation.
        """
        raise NotImplementedError

    @abstractmethod
    def _to_dataset():
        """
        Subclasses must implement this function. `to_dataset()` should
        collect the data from propagate into a dataset which it assigns
        to self.ds.
        """
        raise NotImplementedError

    @abstractmethod
    def _dataset_variables() -> List:
        """
        Should return a list of all of the variables exported in to_dataset(), so
        that has_data() and has_valid_data() check the correct variables.

        This method allows simpler testing.
        """
        raise NotImplementedError


class LimbPropagator(Propagator):
    def __init__(self, platform: Platform, time_offsets: NDArray[timedelta], fov_rotation_matrix: NDArray[float],
                 boresight_angle: float, tangent_altitude_m: float):
        """
        A subclass of Propagator designed for propagating an instrument with
        a limb viewing geometry.

        Parameters
        ----------
        platform : Platform
            An instance of a Skplatform Platform class. The platform should
            have been created by the Instrument class.
        time_offsets : NDArray[timedelta]
            An array of timedeltas, each one indicating the time from TLE
            epoch of an observation.
        fov_rotation_matrix : NDArray[float]
            An array with dimensions [along, across, 3, 3], which has
            rotation matrices for each pixel of the instrument.
        boresight_angle : float
            The angle between 'forward' (the satellite's velocity) and the
            central boresight of the instrument. For a backward looking
            limb instrument this should be 180.0. The angle is in degrees,
            and must be in the range [-360, 360].
        tangent_altitude_m : float
            The altitude the instrument's central boresight looks at, in
            meters. There is a lower limit of 1000 m.
        """
        super().__init__(platform, time_offsets, fov_rotation_matrix)

        self.tangent_altitude_m = 0.0
        self.boresight_angle = 180.0
        self.set_tangent_altitude(tangent_altitude_m)
        self.set_boresight_angle(boresight_angle)

    def set_tangent_altitude(self, value: float):
        """
        Sets the tangent altitude of the propagator. Normally this is only
        called during object creation. Values below 1000 m will raise
        a warning.
        """
        try:
            altitude = float(value)
        except (ValueError, TypeError) as ex:
            logging.error('Trying to set LimbPropagator tangent altitude to a value that could not be converted to a float.')

        if altitude < 1000.0:
            logging.warning(f'Setting LimbPropagator tangent altitude to a value less than 1000.0, {value}, '
                            f'Is the value in kilometers and not meters?')
        self.tangent_altitude_m = altitude

    def set_boresight_angle(self, value: float):
        """
        Sets the boresight angle, in degrees, of the propagator. Normally,
        this is only called during object creation. Values outside the
        range [-360, 360] will raise a warning.
        """
        try:
            angle = float(value)
        except (ValueError, TypeError) as ex:
            logging.error('Trying to set LimbPropagator boresight angle to a value that could not be converted to a float.')

        if (angle < -360.0) or (angle > 360.0):
            logging.warning(f'Setting LimbPropagator boresight angle to a value outside [-360, 360]: {angle}. ')
        self.boresight_angle = angle

    # @njit(fastmath=True, cache=True)
    # this doesn't have time dependence
    def rotate_icf(self) -> NDArray:
        '''
        Applies the rotations for the lines-of-sight to the boresight in the ICF (Instrument Control Frame). In
        the ICF, the boresight is always the x axis, with +x along outward from the instrument.

        Returns
        -------
        NDArray [along_pixels, across_pixels, 3]
            An array of vectors, still in the ICF frame, that represent
            the LOS of each pixel in the detector.
        '''
        # these are to tell njit it can assume these sizes when optimizing
        assert self.R.shape[2] == 3
        assert self.R.shape[3] == 3

        icf_look_vectors = np.zeros((self.num_measurements, self.num_along, self.num_across, 3))
        icf_boresight = np.array([1.0, 0.0, 0.0])
        for time in range(self.num_measurements):
            for along in range(self.num_along):
                for across in range(self.num_across):
                     lv = self.R[along, across] @ icf_boresight
                     icf_look_vectors[time, along, across, :] = lv / np.linalg.norm(lv)

        return icf_look_vectors

    # approximately 10x faster than einsum but only changes milliseconds to microseconds
    # @njit(fastmath=True, cache=True)
    def convert_to_ecef(self, icf_to_ecef_R: NDArray, icf_look_vectors: NDArray) -> NDArray:
        '''
        Converts lines-of-sight in the ICF frame to the ECEF frame.

        Parameters
        ----------
        icf_to_ecef_R : NDArray [num_observations, 3, 3]
            _description_ The 3x3 matrices that transform from the ICF frame to the ECEF frame. These are
            normally calculated by and retrieved from self.platform.icf_to_ecef_rotation_matrices. Each
            observation will have it's own matrix since the satellite orientation and position change
            between observations. All LOS in an observation share the same matrix.

        icf_look_vectors : NDArray [along_pixels, across_pixels, 3, 3]
            _description_ The 3x3 matrices that transform the central line-of-sight (or boresight) to the
            lines-of-sight for each pixel on the detector. These are expected to be provided by
            self.rotate_icf(). There is no time dependence.

        Returns
        -------
        NDArray [num_observations, along_pixels, across_pixels, 3]
            _description_ A multi-dimensional array with the LOS vectors of each pixel on the detector for each
            observation, in the ECEF frame.
        '''
        # these are to tell njit it can assume these sizes when optimizing
        assert icf_look_vectors.shape[3] == 3
        assert icf_to_ecef_R.shape[1] == 3
        assert icf_to_ecef_R.shape[2] == 3

        ecef_look_vectors = np.zeros((self.num_measurements, self.num_along, self.num_across, 3))
        for time in range(self.num_measurements):
            i2e_R = icf_to_ecef_R[time]
            for along in range(self.num_along):
                for across in range(self.num_across):
                    lv = i2e_R @ icf_look_vectors[time, along, across]
                    ecef_look_vectors[time, along, across, :] = lv / np.linalg.norm(lv)

        return ecef_look_vectors

    def propagate(self):
        """
        Propagates the orbit, and calculates the measurement locations,
        ie. the tangent points, for every observation. The values are
        stored in the class variables.
        """
        # dumb hack to shut numpy up about not being able to represent UTC in datetime64
        measurement_times = np.array([t.replace(tzinfo=None) for t in self.measurement_times])
        # sets the local look vector to point at the correct altitude and angle relative to the satellite velocity
        self.platform.add_measurement_set(measurement_times, ('from_platform',),
                                          ('tangent_from_orbitplane', 'standard',
                                           (self.tangent_altitude_m, self.boresight_angle, 0)))
        optical_measurements = self.platform.make_optical_geometry()  # this clears the arrays.... and resets ICF rotation? b/c pointing is from platform

        # can't rotate the actual look vectors in m.look_vector because they are already in ECEF
        # and rotate about that axis not local up

        # # rotate the local look vector in ICF, local look is always [1,0,0] in ICF
        icf_look_vectors = self.rotate_icf()  # -> [time, angle, xyz]

        # get the rotation matrix, this translates the local look ([1,0,0]) from ICF to ECEF
        icf_to_ecef_R = np.array([m.R for m in self.platform.icf_to_ecef_rotation_matrices])  # -> [time, xyz, xyz']

        ecef_look_vectors = self.convert_to_ecef(icf_to_ecef_R, icf_look_vectors)

        sat_positions = np.ones((self.num_measurements, self.num_along, self.num_across, 3))*np.nan
        for time in range(self.num_measurements):
            pos = optical_measurements[time].observer
            for along in range(self.num_along):
                for across in range(self.num_across):
                    sat_positions[time, along, across, :] = pos

        geo = Geodetic()
        tp_lla = geo.llh_from_xyz(geo.xyz_tangent_point_location(sat_positions, ecef_look_vectors))

        latitude = tp_lla[:, :, :, 0]
        longitude = tp_lla[:, :, :, 1]
        altitude = tp_lla[:, :, :, 2]

        sat_positions = np.array([m.observer for m in optical_measurements])
        observer_lla = geo.llh_from_xyz(sat_positions)  # [time, lla]

        observer_latitude = observer_lla[:, 0]
        observer_longitude = observer_lla[:, 1]
        observer_altitude = observer_lla[:, 2]

        wrapped_longitudes = longitude > 180
        longitude[wrapped_longitudes] = longitude[wrapped_longitudes] - 360.0

        wrapped_longitudes = observer_longitude > 180
        observer_longitude[wrapped_longitudes] = observer_longitude[wrapped_longitudes] - 360.0

        time = [np.datetime64(t.replace(tzinfo=None)) for t in self.measurement_times]
        orbit = np.array([self.platform.platform_locator.orbit_number(t)[0] for t in self.measurement_times])

        self._to_dataset(latitude, longitude, altitude, observer_latitude, observer_longitude, observer_altitude, time, orbit)

    def _dataset_variables(self):
        """
        Returns a list of the names of the variables in the dataset.
        """
        # a list of the variables that should be in the dataset.
        return ['latitude', 'longitude', 'altitude', 'time',
                'observer_latitude', 'observer_longitude', 'observer_altitude', 'orbit']

    def _to_dataset(self, latitude: NDArray, longitude: NDArray, altitude: NDArray,
                    observer_latitude: NDArray, observer_longitude: NDArray, observer_altitude: NDArray,
                    time: NDArray, orbit: NDArray):
        """
        Collects all of the data from a simulated orbit and converts it to
        a dataset, which is assigned to self.ds.
        """
        ds = xr.Dataset({
            'latitude': (['time', 'along', 'across'], latitude),
            'longitude': (['time', 'along', 'across'], longitude),
            'altitude': (['time', 'along', 'across'], altitude),
            'time': np.array(time, dtype='datetime64[ns]'),  # convert to ns to silence warning

            'observer_latitude': (['time'], observer_latitude),
            'observer_longitude': (['time'], observer_longitude),
            'observer_altitude': (['time'], observer_altitude),

            'orbit': (['time'], orbit)
        })

        # assign units
        ds.attrs['start_time'] = ds.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        ds.attrs['end_time'] = ds.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        ds.latitude.attrs['units'] = 'degrees_north'
        ds.longitude.attrs['units'] = 'degrees_east'
        ds.altitude.attrs['units'] = 'meters'
        ds.observer_altitude.attrs['units'] = 'meters'
        ds.observer_latitude.attrs['units'] = 'degrees_north'
        ds.observer_longitude.attrs['units'] = 'degrees_east'
        ds.orbit.attrs['units'] = 'orbit number, counted from start of ephemeris'

        self.ds = ds


class NadirPropagator(Propagator):
    def __init__(self, platform: Platform, time_offsets: NDArray[timedelta], fov_rotation_matrix: NDArray):
        """
        A subclass of Propagator designed for propagating an instrument with
        a nadir viewing geometry.

        Parameters
        ----------
        platform : Platform
            An instance of a Skplatform Platform class. The platform should
            have been created by the Instrument class.
        time_offsets : NDArray[timedelta]
            An array of timedeltas, each one indicating the time from TLE
            epoch of an observation.
        fov_rotation_matrix : NDArray[float]
            An array with dimensions [along, across, 3, 3], which has
            rotation matrices for each pixel of the instrument.
        """
        super().__init__(platform, time_offsets, fov_rotation_matrix)

    @staticmethod
    def construct_nadir_rotation_matrix(local_up, look_vector) -> NDArray:
        """
        Constructs a rotation matrix that rotates the central boresight
        from looking along the satellite velocity, the default orentation,
        to looking straight down.

        z| /y                          /y   x=[1,0,0] -> x=[ 0, 0, 1]
         |/__ x  rotate 90 around y:  /__z  y=[0,1,0] -> y=[ 0, 1, 0]
                                      |     z=[0,0,1] -> z=[-1, 0, 0]
                                      |x

        It should be possible to do this using `skplatform`.
        """
        down_axis = -local_up

        across_axis = np.cross(look_vector, down_axis)
        across_axis /= np.linalg.norm(across_axis)

        along_axis = np.cross(down_axis, across_axis)
        along_axis /= np.linalg.norm(along_axis)

        nadir_mat = np.vstack([down_axis, across_axis, along_axis])
        # transpose
        return nadir_mat.T

    def construct_nadir_look_vectors(self, optical_measurements: Any):
        """
        Constructs nadir looking look vectors. The central boresight should
        be looking straight down, with other lines of sight looking at the
        appropriate angles.

        Uses the look vector calculated by Skplatform but not the local_up
        vector. the local look vector is only uses to define the xz plane.

        It should be possible to do this using only `skplatform`.

        Parameters
        ----------
        optical_measurements : Any
            The results of calling platform.make_optical_geometry(). It
            should contain the the local look vectors calculated by
            `skplatform`, which are used to define the xz plane of the
            ICF in the ECEF frame.
        """
        ecef_look_vectors = np.zeros((self.num_measurements, self.num_along, self.num_across, 3))
        # # This does not work; these vectors do not point in the up we want.
        # # We expect the satellite to be oriented so that the nadir is perpendicular
        # # to the surface tangent but it is not... this might be from using the 'limb'
        # # roll control setting. The standard does not work correctly either.
        # local_up_vectors = [m.local_up for m in optical_measurements]
        #                                                           z| /y
        # I suspect the satellite has some rotation about y:         |/__ x
        # so that the central boresight points at the target.
        sat_geo = Geodetic()

        local_up_vectors = []
        for m in optical_measurements:
            west, south, up = sat_geo.xyz_west_south_up(xyz=m.observer)
            local_up_vectors.append(up)

        boresight_vectors = [m.look_vector for m in optical_measurements]
        for t in range(self.num_measurements):
            nadir_mat = self.construct_nadir_rotation_matrix(local_up_vectors[t], boresight_vectors[t])
            for along in range(self.num_along):
                for across in range(self.num_across):
                    rv = nadir_mat @ self.R[along, across, :, :] @ np.array([1.0, 0.0, 0.0])
                    ecef_look_vectors[t, along, across, :] = rv / np.linalg.norm(rv)
        return ecef_look_vectors

    def find_zero_altitude_intercept(self, measurements, look_vectors: NDArray, precision: float = 1.0e-6):
        """
        Finds the approximate location of the zeros altitude intercept. The
        final altitude should be smaller than precision. The functions
        works by projecting the look vector from the satellite position.
        The initial length of the projection is the altitude of the
        satellite, which should be the correct value for the central LOS.
        The height at the end of the projection is then compared to the
        precision, if the height is less then the altitude is returned,
        otherwise the height is added to the length of the projection
        and the end point of the projection is recalculated.

        Parameters
        ----------
        measurements : Any
            The results of calling platform.make_optical_geometry(). Must
            contain the observer xyz position, locating the observer in
            the ECEF frame.
        look_vectors : NDArray [time, along, across, 3]
            The vectors of each of the LOS of the instrument.
        precision : float, optional
            The maximum height above the geoid that will be returned. All
            altitudes below this are considered close enough to zero. The
            default is 1.0e-6.
        """
        intercepts = np.zeros_like(look_vectors)

        for t in range(self.num_measurements):
            satellite = Geodetic()
            sat_altitude = satellite.llh_from_xyz(measurements[t].observer)[2]
            look_vector = look_vectors[t,along,across]

            for along in range(self.num_along):
                for across in range(self.num_across):
                    guess = Geodetic()
                    da = 0.0
                    guess_llh = guess.llh_from_xyz(measurements[t].observer + look_vector*sat_altitude)
                    # simple iterative solver
                    # get the altitude of the current guess, if not close enough add the altitude to the next guess
                    # this works because at nadir this is the exact value needed, and away from nadir
                    # this will always be shorter due to the angle
                    # there may be an issue due to the geoid? or high angles that don't intersect with the ground
                    iter = 0
                    while guess_llh[2] > precision and iter < 10:
                        da += guess_llh[2]
                        length = sat_altitude + da
                        guess_llh = guess.llh_from_xyz(measurements[t].observer + look_vector*length)
                        iter += 1

                    intercepts[t, along, across, :] = guess_llh

    def propagate(self):
        """
        Propagates the orbit, and calculates the measurement locations,
        ie. the zero altitude intercepts, for every observation. The
        resultant data is  stored in the class variables.
        """
        # velocity is only available for the last time simulated
        # so use the tangent looking at the satellite's height
        # and 0 degrees so it is facing the satellite's direction
        sat_altitude_m = self.platform.platform_locator.lat_lon_height[2] - 1_000.0 # exact produces errors
        sat_altitude_m = 20_000.0

        # dumb hack to shut numpy up about not being able to represent UTC in datetime64
        measurement_times = np.array([t.replace(tzinfo=None) for t in self.measurement_times])
        # use of the satellite altitude slightly, ~1e-6, changes the look vectors and therefore the locations
        self.platform.add_measurement_set(measurement_times,
                                          ('from_platform',),  # position
                                          ('tangent_from_orbitplane', 'standard', (sat_altitude_m, 0.0, 0.0)))  # orientation

        optical_measurements = self.platform.make_optical_geometry()  # this clears the arrays.... and resets ICF rotation? b/c pointing is from platform

        ecef_look_vectors = self.construct_nadir_look_vectors(optical_measurements)

        tp_lla = np.zeros_like(ecef_look_vectors)

        for t in range(self.num_measurements):
            satellite = Geodetic()
            sat_altitude = satellite.llh_from_xyz(optical_measurements[t].observer)[2]

            for along in range(self.num_along):
                for across in range(self.num_across):
                    guess = Geodetic()
                    da = 0.0
                    guess_llh = guess.llh_from_xyz(optical_measurements[t].observer + ecef_look_vectors[t,along,across]*sat_altitude)
                    while guess_llh[2] > 1000.0:
                        da += guess_llh[2]
                        length = sat_altitude + da
                        guess_llh = guess.llh_from_xyz(optical_measurements[t].observer + ecef_look_vectors[t,along, across]*length)

                    tp_lla[t, along, across, :] = guess_llh

        latitude = tp_lla[:, :, :, 0]
        longitude = tp_lla[:, :, :, 1]

        sat_positions = np.array([m.observer for m in optical_measurements])
        geo = Geodetic()
        observer_lla = geo.llh_from_xyz(sat_positions)  # [time, lla]

        observer_latitude = observer_lla[:, 0]
        observer_longitude = observer_lla[:, 1]
        observer_altitude = observer_lla[:, 2]

        wrapped_longitudes = longitude > 180.0
        longitude[wrapped_longitudes] = longitude[wrapped_longitudes] - 360.0

        wrapped_longitudes = observer_longitude > 180.0
        observer_longitude[wrapped_longitudes] = observer_longitude[wrapped_longitudes] - 360.0

        time = [np.datetime64(t.replace(tzinfo=None)) for t in self.measurement_times]
        orbit = np.array([self.platform.platform_locator.orbit_number(t)[0] for t in self.measurement_times])

        self._to_dataset(latitude, longitude, observer_latitude, observer_longitude, observer_altitude, time, orbit)

    def _dataset_variables(self) -> List[str]:
        """
        Returns a list of the names of the variables in the dataset.
        """
        return ['latitude', 'longitude', 'time', 'observer_latitude', 'observer_longitude', 'observer_altitude', 'orbit']

    def _to_dataset(self, latitude: NDArray, longitude: NDArray, observer_latitude: NDArray, observer_longitude: NDArray,
                    observer_altitude: NDArray, time: NDArray, orbit: NDArray):
        """
        Collects all of the data from a simulated orbit and converts it to
        a dataset, which is assigned to self.ds.
        """

        ds = xr.Dataset({
            'latitude': (['time', 'along', 'across'], latitude),
            'longitude': (['time', 'along', 'across'], longitude),
            'time': np.array(time, dtype='datetime64[ns]'),  # convert to ns to silence warning

            'observer_latitude': (['time'], observer_latitude),
            'observer_longitude': (['time'], observer_longitude),
            'observer_altitude': (['time'], observer_altitude),

            'orbit': (['time'], orbit)
        })

        # assign units
        ds.attrs['start_time'] = ds.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        ds.attrs['end_time'] = ds.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        ds.latitude.attrs['units'] = 'degrees_north'
        ds.longitude.attrs['units'] = 'degrees_east'
        ds.observer_altitude.attrs['units'] = 'meters'
        ds.observer_latitude.attrs['units'] = 'degrees_north'
        ds.observer_longitude.attrs['units'] = 'degrees_east'
        ds.orbit.attrs['units'] = 'orbit number, counted from start of ephemeris'

        self.ds = ds

class SimPropagator(Propagator):
    def __init__(self, platform: Platform, time_offset: timedelta, fov_rotation_matrix: NDArray,
                 boresight_angle: float, tangent_altitude_m: float):
        """
        Propagator subclass for simulating a limb viewing orbit for the
        `hawcsimulator`. Only a single LOS, the central boresight, is
        calculated, since multiple lines-of-sight are handled by
        `hawcsimulator`.

        Parameters
        ----------
        platform : Platform
            An instance of a Skplatform Platform class. The platform should
            have been created by the Instrument class.
        time_offsets : NDArray[timedelta]
            An array of timedeltas, each one indicating the time from TLE
            epoch of an observation.
        fov_rotation_matrix : NDArray[float]
            An array with dimensions [along, across, 3, 3], which has
            rotation matrices for each pixel of the instrument.
        boresight_angle : float
            The angle between 'forward' (the satellite's velocity) and the
            central boresight of the instrument. For a backward looking
            limb instrument this should be 180.0. The angle is in degrees,
            and must be in the range [-360, 360].
        tangent_altitude_m : float
            The altitude the instrument's central boresight looks at, in
            meters. There is a lower limit of 1000 m.
        """
        super().__init__(platform, time_offset, fov_rotation_matrix)

        self.tangent_altitude_m = 0.0
        self.boresight_angle = 180.0
        self.set_tangent_altitude(tangent_altitude_m)
        self.set_boresight_angle(boresight_angle)

    def set_tangent_altitude(self, value: float):
        """
        Sets the tangent altitude of the progagator. Normally this is only
        called during object creation. Values below 1000 m will raise
        a warning.
        """
        try:
            altitude = float(value)
        except (ValueError, TypeError) as ex:
            logging.error('Trying to set SimPropagator tangent altitude to a value that could not be converted to a float.')

        if altitude < 1000.0:
            logging.warning(f'Setting SimPropagator tangent altitude to a value less than 1000.0, {value}, '
                            f'Is the value in kilometers and not meters?')
        self.tangent_altitude_m = altitude

    def set_boresight_angle(self, value: float):
        """
        Sets the boresight angle, in degrees, of the propagator. Normally,
        this is only called during object creation. Values outside the
        range [-360, 360] will raise a warning.
        """
        try:
            angle = float(value)
        except (ValueError, TypeError) as ex:
            logging.error('Trying to set SimPropagator boresight angle to a value that could not be converted to a float.')

        if (angle < -360.0) or (angle > 360.0):
            logging.warning(f'Setting SimPropagator boresight angle to a value outside [-360, 360]: {angle}. ')
        self.boresight_angle = angle

    def propagate(self):
        """
        Propagates the orbit, and calculates the measurement locations,
        ie. the tangent points, for every observation. The values are
        stored in the class variables.
        """
        # dumb hack to shut numpy up about not being able to represent UTC in datetime64
        # sets the local look vector to point at the correct altitude and angle relative to the satellite velocity
        self.platform.add_measurement_set(self.measurement_times[0].replace(tzinfo=None),
                                        ('from_platform',), ('tangent_from_orbitplane', 'standard',
                                                            (self.tangent_altitude_m, self.boresight_angle, 0)))
        optical_measurement = self.platform.make_optical_geometry()  # this clears the arrays.... and resets ICF rotation? b/c pointing is from platform

        # we don't need to rotate the boresight, there is only one look vector
        icf_boresight = np.array([1.0, 0.0, 0.0])
        # get the rotation matrix, this translates the local look ([1,0,0]) from ICF to ECEF
        icf_to_ecef_R = self.platform.icf_to_ecef_rotation_matrices[0].R
        # convert from ICF to ECEF
        lv = icf_to_ecef_R @ icf_boresight
        ecef_look_vector = lv / np.linalg.norm(lv)

        sat_position = optical_measurement[0].observer

        geo = Geodetic()
        tp_xyz = geo.xyz_tangent_point_location(sat_position, ecef_look_vector)
        tp_lla = geo.llh_from_xyz(tp_xyz)

        time = np.datetime64(self.measurement_times[0].replace(tzinfo=None))
        latitude = tp_lla[0]
        longitude = tp_lla[1]
        altitude = tp_lla[2]

        sat_position = np.array(optical_measurement[0].observer)
        observer_lla = geo.llh_from_xyz(sat_position)  # [time, lla]

        observer_altitude = observer_lla[2]

        if isinstance(longitude, (float, np.float64)):
            if longitude >= 180.0:
                longitude -= 360.0
        else:
            wrapped_longitudes = longitude > 180
            longitude[wrapped_longitudes] = longitude[wrapped_longitudes] - 360.0

        azimuth_angle = self._calculate_north_azimuth_angle(tp_xyz, ecef_look_vector)

        self._to_dataset(latitude, longitude, altitude, observer_altitude, azimuth_angle, time)

    def _calculate_north_azimuth_angle(self, xyz_position: NDArray, look_vector: NDArray ):
        """
        Calculates the angle between a vector and local North at
        the location given by `xyz_position`. The vector is _not_
        projected onto the local tangent plane. At the tangent
        point this does not matter, as the tangent point is when
        the look vector is in the tangent plane, by definition.

        Parameters
        ----------
        xyz_position : NDarray[3]
            The position to calculate North at, in ECEF.
        look_vector : NDArray[3]
            The look vector in ECEF coordinates.
        """
        geo = Geodetic()
        north, east, down = geo.xyz_north_east_down(xyz_position)

        # these should be normed, but play it safe
        unit_north = north / np.linalg.norm(north)
        # Shouldn't this be projected onto the horizontal plane at the tangent point?
        # Generally, yes, but the look vector IS the tangent plane at the tangent point
        unit_lv = look_vector / np.linalg.norm(look_vector)

        # I found two methods for determining the signed angle between two vectors
        # this is the more stable. See tests/angles.py
        cross = np.cross(unit_north, unit_lv)
        up = -down.squeeze()
        unit_up = up / np.linalg.norm(up)
        angle = np.arctan2(np.dot(cross, unit_up), np.dot(unit_lv, unit_north))
        azimuth_angle = np.rad2deg(angle)
        # match the range in hawcsimulator
        if azimuth_angle < 0.0:
            azimuth_angle = 360.0 + azimuth_angle
        return azimuth_angle

    def _dataset_variables(self):
        """
        Returns a list of the names of the variables in the dataset.
        """
        return ['latitude', 'longitude', 'altitude', 'azimuth_angle', 'time', 'observer_altitude']

    def _to_dataset(self, latitude: NDArray, longitude: NDArray, altitude: NDArray, observer_altitude: NDArray, azimuth_angle: NDArray, time: NDArray):
        """
        Collects all of the data from a simulated orbit and converts it to
        a dataset, which is assigned to self.ds.
        """

        ds = xr.Dataset({
            'latitude': latitude,
            'longitude': longitude,
            'altitude': altitude,
            'azimuth_angle': azimuth_angle,
            'time': np.array(time, dtype='datetime64[ns]'),  # convert to ns to silence warning

            'observer_altitude': observer_altitude,
        })

        # assign units
        ds.attrs['date_time'] = ds.time.dt.strftime('%Y-%m-%d %H:%M:%S').item()
        ds.latitude.attrs['units'] = 'degrees_north'
        ds.longitude.attrs['units'] = 'degrees_east'
        ds.altitude.attrs['units'] = 'meters'
        ds.observer_altitude.attrs['units'] = 'meters'
        ds.azimuth_angle.attrs['units'] = 'degrees'

        self.ds = ds
