..  _positioning_technique:

*********************
Positioning Technique
*********************
*Positioning Techniques* are used to specify the position of a platform at a set of measurement times. For example, one technique may specify the position
of a platform using a three element array of geodetic coordinates, :math:`(latitude, longitude, height)`,  while another
technique may specify the same location using a three element array of geocentric coordinates, :math:`(x,y,z)`. The
user chooses which technique they wish to use, see table below.

The following positioning techniques are provided with the ``skplatform`` package,

=============================== ========================================================================================================
Position Technique              Description
=============================== ========================================================================================================
:ref:`xyz`                      The platform is placed at the geocentric ECEF point given by (x,y,z)
:ref:`llh`                      The platform is placed at the geodetic ECEF location given by (latitude, longitude, height)
:ref:`looking_at_llh`           The platform is placed so it is looking at a given location with a given geographic bearing.
:ref:`from_platform`            The platform is placed at the position of a real instrument or platform simulator using :class:`~.PlatformLocation`.
=============================== ========================================================================================================

The last entry, :ref:`from_platform`, allows users to fetch the platform position at a given universal time from an attached
:class:`~.PlatformLocation` object. This technique is used to model various satellite trajectories,
both real and simulated, and is useful when users wish to simulate measurements as a satellite, aircraft or balloon, moves
along its trajectory.

Coding up the positioning technique consists of the following steps

#. Determine the number of measurements. This is driven by the number of different platform positions and orientations
#. Create an array of universal times with one time for each measurement. The times can be the same or they can be different
#. Create an array of platform position arrays. One array entry is created for each measurement.
#. Create a tuple whose first element identifies the position technique and the second element is the array of positions.
#. Pass the array of UTC times and the tuple of position technique into :meth:`~.platform.add_measurement_set` along with orientation info described below.

------------------------------------------------------------------------------------------------------------------------

..  _xyz:

xyz
===
The ``xyz`` positioning technique places the platform at the *N* specified geocentric locations for each measurement entry in
the positioning tuple.

..  function:: Positioning Tuple ('xyz', position)

    :param position[N,3]:
        is N entries of three element arrays,. It can be specified as a sequence of three element arrays, ``[array[3], array[3], ...]``,
        or as an numpy array of shape(N,3). It is internally coerced into a numpy array of dimension (N,3). All values are
        expressed in meters in the :ref:`ecef`.
    :param float [0]:
        **x** component of the platform ECEF position in meters.
    :param float [1]:
        **y** component of the platform ECEF position in meters.
    :param float [2]:
        **z** component of the platform ECEFposition in meters.


Example::

    from skplatform import Platform

    def make_geometry():
        platform = Platform()
        utc = ['2020-09-24T12:15:36.123456', '2020-09-24T12:15:37.456123', '2020-09-24T12:15:38.654321', '2020-09-24T12:15:39.654321']
        positioning_values = [(7223456.0, 1023456.0, 1423456.0), (7523456.0, 923456.0, 1523456.0), (7223456.0, 1023456.0, 1423456.0), (7523456.0, 923456.0, 1523456.0)]
        pointing_values = [(35000, 10), (27000, 5), (24000, 0, 0), (21000, -5, 0)]

        platform.add_measurement_set(utc, ('xyz', positioning_values), ('tangent_altitude', 'limb', pointing_values))
        obspolicy = platform.make_observation_policy()

------------------------------------------------------------------------------------------------------------------------

..  _llh:

llh
===
The ``llh`` positioning technique places the platform at the geodetic location given by latitude, longitude and height in each
entry of the positioning tuple.

..  function:: Positioning Tuple ('llh', values[N,3])


    :param values[N,3]:
        is N entries of three element arrays,. It can be specified as a sequence of three element arrays, ``[array[3], array[3], ...]``,
        or as an numpy array of shape(N,3). It is internally coerced into a numpy array of dimension (N,3). All values are
        expressed in meters in the :ref:`ecef`.
    :param float [0]:
        **Latitude**. The geodetic latitude of the platform position in degrees.
    :param float [1]:
        **Longitude**. The geodetic longitude of the platform position in degrees, positive East.
    :param float [2]:
        **Height**. The height of the platform above sea-level in meters.

Example::

    from skplatform import Platform

    def make_geometry():
        platform = Platform()
        utc = ['2020-09-24T12:15:36.123456', '2020-09-24T12:15:37.456123', '2020-09-24T12:15:38.654321', '2020-09-24T12:15:39.654321']
        positioning_values = [(52, -107, 600000), (53, -107, 600001), (54, -107, 600002), (55,-107, 600003)]
        pointing_values = [(35000, 10, 0), (27000, 5, 0), (24000, 0, 0), (21000, -5, 0)]

        platform.add_measurement_set(utc, ('llh', positioning_values), ('tangent_altitude', 'limb', pointing_values))
        platform.make_complete_measurement_set()

------------------------------------------------------------------------------------------------------------------------

..  _looking_at_llh:

looking_at_llh
==============
The ``looking_at_llh`` positioning technique places the platform at the location in space that looks at the target location.
along the given bearing. Each entry in the positioning tuple requires 5 numbers, described below. It is useful when
placing a platform to look at a known target. The technique is only accurate to ~10 centimeters due to approximations
made in the algorithm.

..  function:: Positioning Tuple ('looking_at_llh', values[N,5])

    :param values[N,5]:
        N entries of five element arrays,  It can be a sequence of 5 element entries, ``[array[5], array[5], ...]``,
        or a two dimensional array of shape [N,5]. It is internally coerced into a numpy array of dimension (N,5).
    :param float [0]:
        **Tangent Latitude**. The latitude of the tangent location in degrees
    :param float [1]:
        **Tangent Longitude**. The longitude of the tangent location in degrees
    :param float [2]:
        **Tangent Altitude**. The height of the tangent above sea-level in meters.
    :param float [3]:
        **Geographic Bearing**. The geographic bearing in degrees of the target from the observer's location. This is
        calculated at the observer's location. 0 is North, 90 is East, 180 is South and 270 is West.
    :param float [4]:
        **Observer Height**. The height of the observer in meters above sea-level.

**Example**. In the following example a satellite is positioned so it is at 600 km altitude and is looking at a bearing of
45 degrees (NE) towards a tangent point at 10 km altitude above (52N, -107E).  A line of sight is chosen so the instrument
is looking in the same geographic bearing at a tangent altitude at 5 km::

    from skplatform import Platform

    platform = Platform()
    utc = ['2020-09-24T18:00:00.0000000']
    positioning_values = [52.0, -107.0, 10000.0, 45.0, 600000.0]
    pointing_values = [5000.0, 45.0]
    platform.add_measurement_set(utc, ('looking_at_llh', positioning_values), ('tangent_altitude', 'limb', pointing_values))
    obspolicy = platform.make_observation_policy()

------------------------------------------------------------------------------------------------------------------------

..  _from_platform:

from_platform
=============
The ``from_platform`` positioning technique uses an instance of :class:`~.PlatformLocation`, which has been given to the
:class:`~.Platform` instance, to calculate the position of the platform at the time of each measurement. This technique is used for
satellite orbit predictors as well as aircraft, balloon and satellites that provide tables of position and time.
The technique is invalid if no :class:`~.PlatformLocation` object is given to the :class:`~.Platform` instance, which is normally
done in the :class:`~.Platform` constructor.

..  function:: Positioning Tuple ('from_platform',)

    No parameters are required. Any parameters passed in are ignored.

Example::

    from skplatform import Platform

    def make_geometry():
        kepler = SatelliteKepler('2020-09-24T12:00:00.000000', period_from_altitude = 600000.0, inclination_radians= radians(97.0), longitude_of_ascending_node_degrees = 82.0, eccentricity= 0.05)
        platform = Platform(platform_locator=kepler)
        utc  = ['2020-09-24T12:15:36.123456', '2020-09-24T12:15:37.456123', '2020-09-24T12:15:38.654321', '2020-09-24T12:15:39.654321']
        pointing_values = [(35000, 10, 0), (27000, 5, 0), (24000, 0, 0), (21000, -5, 0)]
        platform.add_measurement_set(utc, ('from_platform',), ('tangent_altitude', 'limb', pointing_values))
        obs_policy = platform.make_observation_policy()

..  seealso::
    Class ER2, class Carmen, class GroundSite, class SatelliteLocator.
