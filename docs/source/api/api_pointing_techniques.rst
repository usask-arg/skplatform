
..  _pointing_technique:

*****************************
Platform Pointing Techniques
*****************************
Various pointing techniques have been developed that rotate the entire :ref:`pcf` so the :math:`\hat{x}_{PCF}`
unit vector is pointing in the desired direction. Most techniques are useful for instrument simulation work where the user
wishes to have a simulation look at a given location while other techniques are useful when working with real instruments
where attitude solutions are already provided for the platform. All of these techniques only modify the
platform rotation matrix which allows the user to create multiple instrument look vectors in the
:ref:`icf` and rotate them all, in a consistent fashion, to the final reference frame.

The plaform pointing techniques provide a standard set of core methods and users are free to develop
their own solutions if they wish. We provide techniques to support common geometries encountered in atmospheric remote
sensing, :ref:`rollcontrol_limb`, :ref:`rollcontrol_nadir` and :ref:`rollcontrol_standard`. The various
techniques are summarized in the table below. Note that we recommended the :ref:`rollcontrol` setting to be used with
each technique but this is only a guideline and users may choose to use other :ref:`rollcontrol` values if it meets their needs.
This list of techniques is not extensive by any means and we have built the :class:`OrientationTechniques` class so the new
techniques can be added as needs arise, including techniques developed by users.

=============================== =============================== ============================================================================================
Pointing Technique              Recommended Roll Control        Description
=============================== =============================== ============================================================================================
:ref:`tangent_xyz_look`         :ref:`rollcontrol_limb`         Look towards a tangent point with an (x,y,z) unit vector.
:ref:`tangent_altitude`         :ref:`rollcontrol_limb`         Look towards a tangent point with given height and geographic bearing
:ref:`att_location_xyz`         :ref:`rollcontrol_nadir`        Look in nadir towards location given by x,y,z
:ref:`att_location_llh`         :ref:`rollcontrol_nadir`        Look in nadir towards location at geodetic lat, long, height.
:ref:`att_from_platform`         ---                            Set platform orientation from :class`~.Platform` class. This is used for real instruments
:ref:`att_unit_vectors`          ---                            Set platform orientation with explicit :math:`\hat{x}_{ICF}` and :math:`\hat{z}_{ICF}` unit vectors.
:ref:`att_elev_azi_roll`        :ref:`rollcontrol_standard`     Look in azimuth, elevation and roll direction.
:ref:`yaw_pitch_roll`           :ref:`rollcontrol_standard`     Look in yaw, pitch, roll direction
:ref:`tangent_orbplane`         :ref:`rollcontrol_limb`         Look towards a tangent point with a given height and bearing from the orbit plane.
=============================== =============================== ============================================================================================

**Extreme Cases**

It is difficult to provide sensible analysis for extreme cases. For example, the roll control zero point is undefined in
:ref:`rollcontrol_limb` mode when looking directly downwards. Similar conditions occur when looking horizontally in
:ref:`rollcontrol_nadir` mode. The software does detect these extreme conditions and attempts to do something reasonable
but we strongly recommend the user only use :ref:`rollcontrol_limb` and :ref:`rollcontrol_nadir` roll control values for sensible
and nadir geometries.

The :ref:`rollcontrol_standard` setting has difficulties when looking straight up or down if used with
techniques intended for :ref:`rollcontrol_nadir` or :ref:`rollcontrol_limb` geometries as the required azimuth of the instrument :math:`\hat{y}_{ICF}` becomes
undefined. The problem does not exist for the techniques specifically recommended for :ref:`rollcontrol_standard` roll control
as we ensure the users explicitly provides the azimuth information.

------------------------------------------------------------------------------------------------------------------------

..  _tangent_xyz_look:

tangent_xyz_look
================
Configures the platform so the instrument boresight, :math:`\hat{x}_{ICF}`, points in the specified look direction. The
tangent point of the look vector is used as the target location for determining :ref:`rollcontrol` in either
:ref:`rollcontrol_limb` or :ref:`rollcontrol_nadir` modes.

..  function:: ( 'tangent_xyz_look',  roll_control,  parameters )

    ``parameters`` is an N element array of 3 or 4 element arrays, ``[array[3], array[3], ...]`` or ``[array[4], array[4], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,3) or (N,4). All look vectors are expressed in the :ref:`ecef`

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements. Most users will use `limb`.
    :param float [0]:
        **x**. The x component of the look direction unit vector
    :param float [1]:
        **y**.  The y component of the look direction unit vector
    :param float [2]:
        **z**.  The z component of the look direction unit vector
    :param float [3]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*

Example::

    def configure_look( platform: Platform ):
        utc  = ['2020-09-24T12:15:36.123456', '2020-09-24T12:15:37.456123', '2020-09-24 12:15:38.654321']
        pos  = [(52, -107, 600000), (52, -107, 600000), (54, -107, 600002)]
        look = [(0.58311235, -0.43069497, -0.68882642,  0.0), ( 0.58320668, -0.43220012, -0.68780304, 0.0), (0.5833907, -0.43522548, -0.68573615, 0.0)]

        platform.add_measurement_set('limb', utc, ('llh', pos), ('tangent_xyz_look', 'limb', look))
        obspolicy = platform.make_observation_policy()


..  note::
    This technique is intended to be used for tangent altitudes that are below the observer's position and above 5 km
    below sea-level. Look vectors outside this range are discarded.

------------------------------------------------------------------------------------------------------------------------

..  _tangent_altitude:

tangent_altitude
================
The platform looks at the tangent height, geographic bearing and roll specifie in the parameters. This technique should be used
to get the instrument boresight look

..  function:: ( 'tangent_altitude',  roll_control,  parameters )

    ``parameters`` is an N element array of 2 or 3 element arrays, ``[array[2], array[2], ...]`` or ``[array[3], array[3], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,3) or (N,4).

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements. Most users will use `limb`.
    :param float [0]:
        **height**. The height in meters above sea-level of the requested tangent altitude. This value should be less than the observers height
        and greater than 5 km below ground.
    :param float [1]:
        **bearing**.  The geographic, compass bearing in degrees of the tangent direction measured at the observer's location. 0 is North, 90 is East, 180 is South and 270 is West.
    :param float [2]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*


------------------------------------------------------------------------------------------------------------------------

..  _tangent_orbplane:

tangent_from_orbitplane
=======================
The platform looks at the tangent height, bearing and roll specified in the parameters. The bearing is measured from the
orbit plane rather than geographic North. This technique is only available on platforms that have a valid *platform_locator*
object which calculates *velocity* in addition to *position*.

..  function:: ( 'tangent_from_orbitplane',  roll_control,  parameters )

    ``parameters`` is an N element array of 2 or 3 element arrays, ``[array[2], array[2], ...]`` or ``[array[3], array[3], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,3) or (N,4).

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements. Most users will use `limb`.
    :param float [0]:
        **height**. The height in meters above sea-level of the requested tangent altitude. This value should be less than the observers height
        and greater than 5 km below ground.
    :param float [1]:
        **bearing**.  The bearing in degrees of the tangent point from the orbit plane. The forward looking direction,
        *parallel* to the platform velocity is the point of zero bearing. Bearing increases in the same direction as a
        compass bearing, i.e. clockwise when viewed from above: North->East->South->West.
    :param float [2]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*

..  note::
    The *platform_locator* can be set when you create an instance of class :class:`~.Platform`.

------------------------------------------------------------------------------------------------------------------------

..  _att_location_xyz:

location_xyz
============
The platform looks at the given geocentric location (x,y,z) with roll of 'roll' degrees. Typically used for satellite nadir observations

..  function:: ( 'location_xyz',  roll_control,  parameters )

    ``parameters`` is an N element array of 3 or 4 element arrays, ``[array[3], array[3], ...]`` or ``[array[4], array[4], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,3) or (N,4). All position vectors are expressed in meters in the :ref:`ecef`

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements.
    :param float [0]:
        **x**. The x component of the target position vector
    :param float [1]:
        **y**.  The y component of the target position vector
    :param float [2]:
        **z**.  The z component of the target position  vector
    :param float [3]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*


------------------------------------------------------------------------------------------------------------------------

..  _att_location_llh:

location_llh
============
The platform looks at the given geodetic location (lat, lng, height) with roll of 'roll' degrees. Typically used for satellite nadir observations

..  function:: ( 'location_llh',  roll_control,  parameters )

    ``parameters`` is an N element array of 3 or 4 element arrays, ``[array[3], array[3], ...]`` or ``[array[4], array[4], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,3) or (N,4). All position vectors are expressed in meters in the :ref:`ecef`

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements.
    :param float [0]:
        **Latitude**. The geodetic latitude of the target position in degrees.
    :param float [1]:
        **Longitude**. The geodetic longitude of the target position in degrees
    :param float [2]:
        **Height**. The height of the target position above sea-level in meters.
    :param float [3]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*


------------------------------------------------------------------------------------------------------------------------

..  _att_elev_azi_roll:

azi_elev
========
The platform looks in the direction given by azimuth, elevation and roll (applied in that order). Used for ground sites


..  function:: ( 'azi_elev',  roll_control,  parameters )

    ``parameters`` is an N element array of 2 or 3 element arrays, ``[array[2], array[2], ...]`` or ``[array[3], array[3], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,2) or (N,3).

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements.
    :param float [0]:
        **Azimuth**. The geographic azimuth of the instrument boresight in degrees from North. Measured clockwise from North. 0 is North, 90 is East, 180 is South, 270 is West.
    :param float [1]:
        **Elevation**. The elevation in degrees of the instrument boresight from the horizontal plane at the observer's location. Positive elevation (0-90) is upwards. Negative elevation (-90 to 0) is downwards.
    :param float [2]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*. Most users will use *standard*.

-----------------------------------------------------------------------------------------------------------------------

..  _yaw_pitch_roll:

yaw_pitch_roll
==============
The platform applies pointing information in the order yaw, pitch, roll. This is useful for aircraft and balloon systems.
Most users will choose to use *standard* :ref:`rollcontrol`.


..  function:: ( 'yaw_pitch_roll',  roll_control,  parameters )

    ``parameters`` is an N element array of 2 or 3 element arrays, ``[array[2], array[2], ...]`` or ``[array[3], array[3], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,2) or (N,3).

    :param str roll_control:
        The :ref:`rollcontrol` value applied to this set of measurements.
    :param float [0]:
        **Yaw**. The geographic bearing of the instrument boresight in degrees from North. Measured clockwise from North. 0 is North, 90 is East, 180 is South, 270 is West.
    :param float [1]:
        **Pitch**. The pitch or elevation elevation in degrees of the instrument boresight from the horizontal plane at the observer's location. Positive elevation (0-90) is upwards. Negative elevation (-90 to 0) is downwards.
    :param float [2]:
        **roll**. Optional [default 0.0]. The roll angle in degrees of the instrument control frame around the boresight,
        :math:`\hat{x}_{ICF}`, from the zero point implied by *roll_control*. Most users will use *standard*.

------------------------------------------------------------------------------------------------------------------------

..  _att_from_platform:

from_platform
=============
The platform orientation is set by the parameters returned by the :class`~.Platform` class at the required times. This is used for real instruments

..  function:: ( 'from_platform',) or ('from_platform')

    No limb control values or parameters are required.


------------------------------------------------------------------------------------------------------------------------

..  _att_unit_vectors:

unit_vectors
============
The platform orientation is explicitly  so the instrument :math:`\hat{x}_{ICF}` and :math:`\hat{z}_{ICF}` are positioned along
the two unit-vectors defined by the 6 parameters. The orientation of :math:`\hat{y}_{ICF}` is given by forming a right-handed
orthogonal system. The *roll_control* value is ignored. Users should be careful as no checks are made to ensure the input
vectors are orthogonal unit vectors.

..  function:: ( 'unit_vectors',  roll_control,  parameters )

    ``parameters`` is an N element array of 6 element arrays, ``[array[6], array[6], ...]``.
    It is anything that can be sensibly coerced into a numpy array of dimension (N,6).

    :param str roll_control:
        This parameter is ignored.
    :param float [0]:
        **Xx**. The x component of :math:`\hat{x}_{ICF}` in the :ref:`ecef`
    :param float [1]:
        **Xy**.  The y component of :math:`\hat{x}_{ICF}` in the :ref:`ecef`
    :param float [2]:
        **Xz**.  The z component of :math:`\hat{x}_{ICF}` in the :ref:`ecef`
    :param float [0]:
        **Zx**. The x component of :math:`\hat{z}_{ICF}` in the :ref:`ecef`
    :param float [1]:
        **Zy**.  The y component of :math:`\hat{z}_{ICF}` in the :ref:`ecef`
    :param float [2]:
        **Zz**.  The z component of :math:`\hat{z}_{ICF}` in the :ref:`ecef`


