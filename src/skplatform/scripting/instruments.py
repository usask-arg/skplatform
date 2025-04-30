from typing import Tuple, List, Optional, Union
import numpy as np
import xarray as xr
import logging
from numpy.typing import NDArray
from datetime import timedelta

from skplatform.platform import Platform
from skplatform.satellite import SatelliteSGP4

from .propagator import Propagator, LimbPropagator, NadirPropagator, SimPropagator
from .rotations import calculate_2d_rotation_matrices



# '1 99999U 99999 A  19213.00000000  .00000000  00000-0  00000-0 0    10\n'
# '2 99999 097.2222 331.8182 0001962 091.0782 359.0176 15.37843164    16\n'

class PropagatorInstrument():
    '''
    The PropagatorInstrument class is the base class for other Instrument classes.
    It is not meant to be used directly. The PropagatorInstrument sub-classes are
    intended to be the user facing interface of the simulator. They are meant to
    manage:
    * the properties of the detector: field-of-view, number of pixels and
    observation frequency
    * the Propagator class associated with the instrument (Limb or Nadir)
    * the simulation time

    The base class has functions to ensure the propagator is configured and updated
    correctly according to user input. Custom subclasses should provide support for
    any custom properties. If the rotation matrices are invalidated by a change
    `self.pixel_grid_rotation_matrix` should be set to None. If the Propagator is
    invalidated `self._propagator_out_of_date` should be set to True. Finally, the
    subclass must implement `self._build_propagator()`. This function should check
    `self.pixel_grid_rotation_matrix` and `self._propagator_out_of_date` and take
    appropriate action if either is invalid.

    The classes `LimbInstrument` and `NadirInstrument` are standard implementations. They
    are provided for custom instruments that only require custom detector parameters.
    '''
    def __init__(self, tle: List[str], along_fov: float, across_fov: float, along_pixels: int, across_pixels: int, frequency: Union[NDArray, float]):
        """
        tle : List[str]
            A two-line element set loaded as two strings. The first line
            should be the first element and the second line is the second
            element.
        along_fov : float
            The field of view in 'forward' direction of the instrument,
            aligned with the satellite ground track by default. The field
            of view is centered, so the angles extend from -along_fov/2 to
            +along_fov/2.
        across_fov : float
            The field of view in 'sideways' direction of the instrument,
            perpendicular with the satellite ground track by default. The
            field of view is centered, so the angles extend from
            -across_fov/2 to +across_fov/2.
        along_pixels : int
            The number of pixels in the along_fov. Each pixel is a angle.
        across_pixels : int
            The number of pixels in the across_fov. Each pixel is a angle.
        frequency : Union[NDArray, float]
            Specifies the interval, in seconds, between measurements. Can be
            either a single float or an array of floats. If a float, every
            observation is frequency seconds apart.

            If an array, the frequency defines a series of observations, an
            observation set. The set is then repeated over the length of the
            simulation. The first element in the set is the interval between
            the first and second observation, and the last element is the
            interval between the last observation in the set, and the first
            observation of the next set.
        """
        self.tle: List[str] = tle
        self.sat = self._setup_satellite()
        self.platform: Platform = Platform(platform_locator=self.sat)
        # use properties?
        # if any of these are set the rotation matrices need to be recalculated
        self._along_fov: float = along_fov
        self._across_fov: float = across_fov
        self._along_pixels: int = along_pixels
        self._across_pixels: int = across_pixels

        self._frequency: float = frequency # Hz

        # self.pixel_grid_rotation_matrix: NDArray[np.float64] = np.array([])
        # delaying calculating these allows the user to change any of the input values
        # without worrying about anything. The only reason to precalculate them is to
        # inspect/edit them
        self.pixel_grid_rotation_matrix: Optional[NDArray[np.float64]] = \
                calculate_2d_rotation_matrices(self.along_fov, self.across_fov, self.along_pixels, self.across_pixels)

        self.start_time: float = 0.0
        self.end_time: float = 60.0

        self.time_offsets: NDArray = np.array([])
        self.propagator: Propagator = None
        self._propagator_out_of_date: bool = True
        # delay building the propagator until propagate is called
        # this lets users editing the instruments change the
        # properties without worrying about anything
        # self._build_propagator()

    @property
    def along_fov(self):
        return self._along_fov

    @along_fov.setter
    def along_fov(self, val: float):
        self._along_fov = val
        self.pixel_grid_rotation_matrix = None  # need to recalculate the grid now
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix

    @property
    def across_fov(self):
        return self._across_fov

    @across_fov.setter
    def across_fov(self, val: float):
        self._across_fov = val
        self.pixel_grid_rotation_matrix = None  # need to recalculate the grid now
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix

    @property
    def along_pixels(self):
        return self._along_pixels

    @along_pixels.setter
    def along_pixels(self, val: float):
        self._along_pixels = val
        self.pixel_grid_rotation_matrix = None  # need to recalculate the grid now
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix

    @property
    def across_pixels(self):
        return self._across_pixels

    @across_pixels.setter
    def across_pixels(self, val: float):
        self._across_pixels = val
        self.pixel_grid_rotation_matrix = None  # need to recalculate the grid now
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix

    # frequency will need adjustment to handle variable frequecies
    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, val: Union[NDArray, float]):
        self._frequency = val

    def __str__(self):
        # define a function so can print these nicely
        txt = (f'FOV Along Track: {self.along_fov}\n'
               f'FOV Across Track: {self.across_fov}\n'
               f'Pixels (LOS) Along Track: {self.along_pixels}\n'
               f'Pixels (LOS) Across Track: {self.across_pixels}\n'
               f'Observation Frequency: {self.frequency} Hz\n'
               f'Simulation Start Time: {self.start_time} s from epoch\n'
               f'Simulation End Time: {self.end_time} s from epoch\n'
               )
        return txt

    def configure_detector(self, along_fov: float, across_fov: float, along_pixels: int, across_pixels: int, frequency: Union[NDArray, float]):
        """
        Reconfigures the instrument properties, except for the tle.

        Parameters
        ----------
        along_fov : float
            The field of view in 'forward' direction of the instrument,
            aligned with the satellite ground track by default. The field
            of view is centered, so the angles extend from -along_fov/2 to
            +along_fov/2.
        across_fov : float
            The field of view in 'sideways' direction of the instrument,
            perpendicular with the satellite ground track by default. The
            field of view is centered, so the angles extend from
            -across_fov/2 to +across_fov/2.
        along_pixels : int
            The number of pixels in the along_fov. Each pixel is a angle.
        across_pixels : int
            The number of pixels in the across_fov. Each pixel is a angle.
        frequency : Union[NDArray, float]
            Specifies the interval, in seconds, between measurements. Can be
            either a single float or an array of floats. If a float, every
            observation is frequency seconds apart.

            If an array, the frequency defines a series of observations, an
            observation set. The set is then repeated over the length of the
            simulation. The first element in the set is the interval between
            the first and second observation, and the last element is the
            interval between the last observation in the set, and the first
            observation of the next set.
        """
        self.along_fov: float = along_fov
        self.across_fov: float = across_fov
        self.along_pixels: int = along_pixels
        self.across_pixels: int = across_pixels
        self.frequency: Union[NDArray, float] = frequency # Hz

        self.pixel_grid_rotation_matrix: NDArray[np.float64] = \
                calculate_2d_rotation_matrices(self.along_fov, self.across_fov, self.along_pixels, self.across_pixels)
        self._propagator_out_of_date: bool = True

    def _setup_satellite(self):
        '''
        Constructs the `skplatform` satellite using the loaded TLE. An error
        will be generated if the TLE is invalid for any reason.

        Returns
        -------
        sat : SatelliteSGP4
            A SGP4 orbit predictor object initialized with the classes TLE.
        '''
        if isinstance(self.tle, list) and (len(self.tle)==2) \
                    and isinstance(self.tle[0], str) and isinstance(self.tle[0], str):
            sat = SatelliteSGP4(twolines=self.tle)
            sat.use_polar_motion = True
            return sat
        else:
            logging.error('TLE for PropagatorInstrument class/subclass is incorrect or incomplete. '
                          'self.tle should be: [line1, line2] where line1 and line2 are the TLE lines as strings')

    def _build_propagator(self):
        # Must be implemented in subclasses
        raise NotImplementedError

    def set_start_time(self, start_time: float):
        '''
        Sets the start time for the orbit simulation in number of seconds
        from the TLE epoch. Start time must be earlier than end time.

        Parameters
        ----------
        start_time : float
            The number of seconds from the TLE epoch to begin orbit
            simulations from. Negative values are permitted.
        '''
        self.start_time = start_time
        self._propagator_out_of_date = True

    def set_end_time(self, end_time: float):
        """
        Sets the end time of the orbit simulation, in seconds from the TLE
        epoch. End time must be later than start time.

        Parameters
        ----------
        end_time : float
            The number of second from the TLE epoch at which to end the
            orbit simulation at. Negative values are permitted.
        """
        self.end_time = end_time
        self._propagator_out_of_date = True

    def set_propagation_time(self, start_time, end_time):
        """
        Sets both the start and end times of the orbit simulation, in
        seconds from the TLE epoch. Either/both values can be negative
        but end time must be chronically later than start time.

        Parameters
        ----------
        start_time : float
            The number of seconds from the TLE epoch to begin orbit
            simulations from. Negative values are permitted.
        end_time : float
            The number of second from the TLE epoch at which to end the
            orbit simulation at. Negative values are permitted.
        """
        # to do enforce start_time < end_time
        self.set_start_time(start_time)
        self.set_end_time(end_time)
        self._propagator_out_of_date = True

    def propagate(self):
        """
        Propagates the orbit over the interval from start time to end
        time. The results are stored in the appropriate object attributes.

        If the propagator was invalidated in any way it will be rebuilt.
        """
        if self._propagator_out_of_date:
            self._build_propagator()
        self.propagator.propagate()

    def propagate_for_time(self, start_time, end_time):
        """
        Sets both the start and end times of the orbit simulation, in
        seconds from the TLE epoch. Either/both values can be negative
        but end time must be chronologically later than start time. After
        setting the times the orbit is propagated.

        Parameters
        ----------
        start_time : float
            The number of seconds from the TLE epoch to begin orbit
            simulations from. Negative values are permitted.
        end_time : float
            The number of second from the TLE epoch at which to end the
            orbit simulation at. Negative values are permitted.
        """
        self.set_propagation_time(start_time, end_time)
        self.propagate()

    def calculate_time_offsets(self) -> NDArray[np.float64]:
        """
        Calculates the time offsets for the observations over the time
        interval specified by start time and end time. The time offsets
        are the number of seconds from the TLE epoch. Frequency as a
        constant or an array of intervals are both supported.

        Returns
        -------
        NDArray
            An array of floats specifying when each observation in the
            simulation will take place, in terms of the number of seconds
            from the TLE epoch.
        """
        if isinstance(self.frequency, float):
            time_offsets = np.array([timedelta(seconds=dt) for dt in np.arange(self.start_time,
                                                                        self.end_time,
                                                                        self.frequency)])
        elif isinstance(self.frequency, (np.ndarray, List)):
            obs_set_length = np.sum(self.frequency)
            chunk_length = self.end_time - self.start_time
            num_obs_sets = int(chunk_length//obs_set_length)
            time_offsets = np.unique(np.cumsum(np.tile(self.frequency, num_obs_sets)))
            time_offsets = np.array([timedelta(seconds=self.start_time + dt) for dt in time_offsets[:-1]])
        else:
            logging.error('Frequency was neither a float for a constant measuremnt interval or a List/NDarray of floats.')
            return
        return time_offsets

    def get_dataset(self) -> Optional[xr.Dataset]:
        """
        Returns the data from propagating the orbit as a xarray Dataset.
        Returns None if there is no data.
        """
        # don't want to call propagate(), as that gives the function
        # side effects that are not obvious.
        if self.propagator is not None:
            # check if has ran the propagator, run it if necessary?
            if self.propagator.ds is not None:
                return self.propagator.ds
        return None

class LimbInstrument(PropagatorInstrument):
    def __init__(self, name: str, tle: List[str], along_fov: float, across_fov: float, along_pixels: int, across_pixels: int,
                 frequency: Union[NDArray, float], look_angle_deg: float, target_altitude_m: float):
        """
        A class for simulating generic Limb instruments.

        name : str
            A name used when saving files generated with this instrument.
        tle : List[str]
            A two-line element set loaded as two strings. The first line
            should be the first element and the second line is the second
            element.
        along_fov : float
            The field of view in 'forward' direction of the instrument,
            aligned with the satellite ground track by default. The field
            of view is centered, so the angles extend from -along_fov/2 to
            +along_fov/2.
        across_fov : float
            The field of view in 'sideways' direction of the instrument,
            perpendicular with the satellite ground track by default. The
            field of view is centered, so the angles extend from
            -across_fov/2 to +across_fov/2.
        along_pixels : int
            The number of pixels in the along_fov. Each pixel is a angle.
        across_pixels : int
            The number of pixels in the across_fov. Each pixel is a angle.
        frequency : Union[NDArray, float]
            Specifies the interval, in seconds, between measurements. Can be
            either a single float or an array of floats. If a float, every
            observation is frequency seconds apart.

            If an array, the frequency defines a series of observations, an
            observation set. The set is then repeated over the length of the
            simulation. The first element in the set is the interval between
            the first and second observation, and the last element is the
            interval between the last observation in the set, and the first
            observation of the next set.
        look_angle_deg : float
            The angle between the satellites's forward direction (its
            velocity) and the instrument's forward direction, in degrees.
        target_altitude_m : float
            The altitude of the instrument's central line-of-sight, or
            boresight, in meters.
        """
        super().__init__(tle, along_fov=along_fov, across_fov=across_fov, along_pixels=along_pixels, across_pixels=across_pixels, frequency=frequency)
        self.name: str = name
        self.look_angle_deg: float = look_angle_deg
        self.target_altitude_m: float = target_altitude_m

    def _build_propagator(self):
        """
        Constructs and configures the propagator using the class attributes.

        Should be called when constructing the class or when the propagator
        is invalid.
        """
        if self.pixel_grid_rotation_matrix is None:
            self.pixel_grid_rotation_matrix = calculate_2d_rotation_matrices(self.along_fov, self.across_fov, self.along_pixels, self.across_pixels)

        self.time_offsets = self.calculate_time_offsets()

        self.propagator = LimbPropagator(self.platform, self.time_offsets, self.pixel_grid_rotation_matrix,
                                         self.look_angle_deg, self.target_altitude_m)

    def __str__(self):
        txt = f'Limb Instrument: {self.name}\n'
        txt += super().__str__()
        txt += (f'Look Angle: {self.look_angle_deg} degrees from forward\n'
                f'Target Altitude: {self.target_altitude_m} m\n')
        return txt

    def configure_limb(self, look_angle_deg: float, target_altitude_m: float, ):
        """
        Configures attributes specific to the limb instruments, which are
        not configured in `configure_detector()`

        Parameters
        ----------
        look_angle_deg : float
            The angle between the satellites's forward direction (its
            velocity) and the instrument's forward direction, in degrees.
        target_altitude_m : float
            The altitude of the instrument's central line-of-sight, or
            boresight, in meters.
        """
        self.look_angle_deg = look_angle_deg
        self.target_altitude_m = target_altitude_m
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix


class AOSSky(LimbInstrument):
    def __init__(self, tle: List[str]):
        """
        A Limb instrument configured with properties of the AOS-Sky
        instrument. Note that the details of the instrument have not
        been finalized and the properties may change between versions.
        """
        super().__init__('AOS-sky', tle, along_fov=0, across_fov=5.0, along_pixels=1, across_pixels=512,
                         frequency=1.0, look_angle_deg=180.0, target_altitude_m=20_000.0)
    def __str__(self):
        txt = super().__str__()
        return txt.replace('Generic Limb Instrument', 'AOS-Sky')


class NadirInstrument(PropagatorInstrument):
    """
    A class for simulating generic Nadir instruments.

    name : str
        A name used when saving files generated with this instrument.
    tle : List[str]
        A two-line element set loaded as two strings. The first line
        should be the first element and the second line is the second
        element.
    along_fov : float
        The field of view in 'forward' direction of the instrument,
        aligned with the satellite ground track by default. The field
        of view is centered, so the angles extend from -along_fov/2 to
        +along_fov/2.
    across_fov : float
        The field of view in 'sideways' direction of the instrument,
        perpendicular with the satellite ground track by default. The
        field of view is centered, so the angles extend from
        -across_fov/2 to +across_fov/2.
    along_pixels : int
        The number of pixels in the along_fov. Each pixel is a angle.
    across_pixels : int
        The number of pixels in the across_fov. Each pixel is a angle.
    frequency : Union[NDArray, float]
        Specifies the interval, in seconds, between measurements. Can be
        either a single float or an array of floats. If a float, every
        observation is frequency seconds apart.

        If an array, the frequency defines a series of observations, an
        observation set. The set is then repeated over the length of the
        simulation. The first element in the set is the interval between
        the first and second observation, and the last element is the
        interval between the last observation in the set, and the first
        observation of the next set.
    """
    # doesn't any special initialization
    def __init__(self, name: str, tle: List[str], along_fov: float, across_fov: float, along_pixels: int, across_pixels: int, frequency: Union[NDArray, float]):
        super().__init__(tle, along_fov=along_fov, across_fov=across_fov, along_pixels=along_pixels, across_pixels=across_pixels, frequency=frequency)
        self.name: str = name

    def _build_propagator(self):
        self.time_offsets = self.calculate_time_offsets()
        self.propagator = NadirPropagator(self.platform, self.time_offsets, self.pixel_grid_rotation_matrix)

    def __str__(self):
        txt = f'Nadir Instrument: {self.name}\n'
        txt += super().__str__()
        return txt


class TICFIRE(NadirInstrument):
    # 38x39 with FOV ±4.76° x ±4.9656°  (along x across)
    def __init__(self, tle: List[str]):
        """
        A Nadir instrument configured with properties of the TICFIRE
        instrument. Note that the details of the instrument have not
        been finalized and the properties may change between versions.
        """
        super().__init__('TICFIRE', tle, along_fov= 4.76*2.0, across_fov=4.9656*2.0,
                        along_pixels=38, across_pixels=39, frequency=10.5)
    def __str__(self):
        txt = super().__str__()
        return txt.replace('Generic Nadir Instrument', 'TICFIRE')


class SimulatorInstrument(PropagatorInstrument):
    """
    A class implementing Limb geometry for use in the `hawcsimulator` package.
    """
    def __init__(self, tle: List[str], look_angle_deg: float = 0.0, target_altitude_m: float = 25_000.0):
        """
        An instrument class for use in the `hawcsimulator` package. It
        implements a limb viewing geometry. Detector properties are
        managed by `hawcsimulator`.

        tle : List[str]
            A two-line element set loaded as two strings. The first line
            should be the first element and the second line is the second
            element.
        look_angle_deg : float
            The angle between the satellites's forward direction (its
            velocity) and the instrument's forward direction, in degrees.
        target_altitude_m : float
            The altitude of the instrument's central line-of-sight, or
            boresight, in meters.
        """
        super().__init__(tle, along_fov=0, across_fov=0, along_pixels=1, across_pixels=1, frequency=1.0)

        self.look_angle_deg: float = look_angle_deg
        self.target_altitude_m: float = target_altitude_m

    def _build_propagator(self):
        if self.pixel_grid_rotation_matrix is None:
            self.pixel_grid_rotation_matrix = calculate_2d_rotation_matrices(self.along_fov, self.across_fov, self.along_pixels, self.across_pixels)

        self.time_offsets = self.calculate_time_offsets()

        self.propagator = SimPropagator(self.platform, self.time_offsets, self.pixel_grid_rotation_matrix,
                                         self.look_angle_deg, self.target_altitude_m)

    def __str__(self):
        txt = 'Generic Limb Instrument\n'
        txt += super().__str__()
        txt += (f'Look Angle: {self.look_angle_deg} degrees from forward\n'
                f'Target Altitude: {self.target_altitude_m} m\n')
        return txt

    def configure_limb(self, look_angle_deg: float, target_altitude_m: float, ):
        """
        Configures attributes specific to the limb instruments, which are
        not configured in `configure_detector()`

        Parameters
        ----------
        look_angle_deg : float
            The angle between the satellites's forward direction (its
            velocity) and the instrument's forward direction, in degrees.
        target_altitude_m : float
            The altitude of the instrument's central line-of-sight, or
            boresight, in meters.
        """
        self.look_angle_deg = look_angle_deg
        self.target_altitude_m = target_altitude_m
        self._propagator_out_of_date = True  # rebuilding the propagator, checks the rotation matrix
