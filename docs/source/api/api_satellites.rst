==================
Satellites
==================
Satellite classes are used to locate the position of a platform for a variety of orbits. They are derived from class :class:`~.PlatformLocation`
which allows them to be used as ``platform locators`` in class :class:`~.Platform` when coupled with the :ref:`from_platform`
positioning technique.

Example::

    from skplatform import Platform

    def make_geometry():
        kepler = SatelliteKepler('2020-09-24T12:00:00.000000', period_from_altitude = 600000.0, inclination_radians= radians(97.0), longitude_of_ascending_node_degrees = 82.0, eccentricity= 0.05)
        platform = Platform(platform_locator=kepler)
        utc  = ['2020-09-24T12:15:36.123456', '2020-09-24T12:15:37.456123', '2020-09-24T12:15:38.654321', '2020-09-24T12:15:39.654321']
        pointing_values = [(35000, 10, 0), (27000, 5, 0), (24000, 0, 0), (21000, -5, 0)]
        platform.add_measurement_set(utc, ('from_platform',), ('tangent_altitude', 'limb', pointing_values))

The satellite classes can also be used as stand-alone classes in which case the user would execute the following steps.

#. Create a specific instance of a satellite.
#. Update the satellite position to a given instant in time, see :meth:`~.PlatformLocation.update_position`
#. Get the satellite position using one of the available methods, see :meth:`~.PlatformLocation.position`, :meth:`~.PlatformLocation.earth_location` or :meth:`~.PlatformLocation.lat_lon_height`
#. Repeat steps 2 and 3.

for example::

    dt = timedelta(minutes=1.0)                                         # step time of our simulation is one minute
    numtimes = 1440                                                     # Simulate for 1 day = 1440 * one minute
    platform_utc = datetime(2019,7,26, hour =20, minute=15, second=00)  # start time of our simulation
    sat = SatelliteSunSync(platform_utc,                                # Create the satellite
                           orbittype='sgp4',
                           period_from_altitude=600000.0,
                           localtime_of_ascending_node_hours=18.25)

    answers = np.zeros([3, numtimes])                                   # make an array to hold the latitude, longitude, altitude of each satellite at each simulation step
    times = np.zeros([numtimes], dtype='datetime64[us]')                # make an array to hold the UTC time of each simulation step.
    for i in range(numtimes):                                           # for each time step
        tnow = platform_utc + i * dt                                    # get the next time step
        times[i] = tnow                                                 # save the time of the step
        sat.update_position(tnow)                                       # Update the satellite position
        answers[:, i] = sat.lat_long_height                             # convert XYZ position and save

Users should be aware the satellite classes may internally implement numerical integration using Runge-Kutta algorithms and
should avoid large time steps (e.g. multiple orbits) as this may result in long execution times in the best case and
inaccuracy in the worst case.

When greater accuracy is required, the *use_polar_motion* flag can be set to True, which corrects for the changes in Earth's rotational axis. 
The first time this is used an internet connection is required to download the data file containing the corrections. The gain in accuracy is 
small and is not expected to be required in normal use cases. All satellite classes have this flag. 


=====================
Satellite Propagators
=====================

..  _satellitekepler_class:

SatelliteKepler
===============
Implements a classic Kepler orbit which generates elliptical orbits with the centre of the Earth located at one focus.
This is quick and fast but keep in mind that it has no gravitational perturbations so it can never truly mimic
sun-synchronous or other classes of orbit that depend up gravitional torques etc.

Example::

    from skplatform.satellite  import SatelliteKepler

    kepler  = SatelliteKepler ( platform_utc, period_from_altitude = 600000.0, inclination_radians= radians(97.0), longitude_of_ascending_node_degrees = 82.0, eccentricity= 0.05)
    pos = kepler.update_position( tnew )

.. autoclass:: skplatform.satellite.SatelliteKepler
    :members:
    :special-members: __init__


..  _satellitesgp4_class:

SatelliteSGP4
=============
A satellite class that uses the SGP4 satellite propagation model implemented in python package `sgp4`.  The classic
usage is to initialize the code with two line elements that you have obtained from `Space Trak <https://www.space-track.org/>`_
or `Celestrak <https://celestrak.com/>`_ or elsewhere. Example using two line elements::

    platform_utc     = datetime(2019,7,26, hour =20, minute=15, second=00)

    line1   = "1 26702U 01007A   19206.65122582  .00000348  00000-0  26798-4 0  9995"
    line2   = "2 26702  97.5720 223.8268 0009640 316.2599  43.7864 15.07871693  7426"
    sgp4    = SatelliteSGP4   ( twolines= [line1,line2] )

The class also provides the option to initialize the SGP4 propagator with the state vector and orbital parameters from a
Kepler orbit. This feature is used to implement some of the other satellite classes (eg sun-sync) in `skretrieval` and provides
a convenient way to use the gravitational perturbations in the SGP4 model without needing a set of two line elements.

.. autoclass:: skplatform.satellite.satellitesgp4.SatelliteSGP4
    :members:
    :special-members: __init__

SatelliteSunSync
================
Create a satellite which is sun-synchronous. The class is intended for modelling/simulation work where a sun-synchronous
orbit is required. The code simulates the sun-synchronous orbit by setting the initial orbit inclination, it is the underlying
orbit propagator's job to calculate the orbit precession. The 'sgp4' orbit propagator will precess the orbit properly at
one degree per day. A simple 'kepler' orbit propagator can also be used but will not precess the orbit, this is satisfactory
for short period simulations over a few hours but should not be used for multi-day or longer studies.

For example, the following code will create a sun-synchronous orbit at a *nominal* altitude of 600 kms awith an ascending node at 18:15LT::

    sunsync = SatelliteSunSync( platform_utc, orbittype='sgp4',
                                period_from_altitude=600000.0,
                                localtime_of_ascending_node_hours=18.25)


Regardless of orbit propagator chosen, the sun synchronous orbit is kick-started with a state vector generated by a Kepler orbit
with the appropriate sun-synchronous inclination.  The state vector is fed to the orbit propagator which then calculates all future
platform_ecef_positions of the satellite.

.. autoclass:: skplatform.satellite.SatelliteSunSync
    :members:
    :special-members: __init__

SatelliteMolniya
================
Create a satellite which follows a Molniya orbit. The class is intended for modelling/simulation work and can use either
a simple Kepler orbit or an SGP4 orbit. For example, the following code will create a molniya orbit::

    molniya   = SatelliteMolniya(  self._utc,
                                   orbittype='sgp4',
                                   longitude_of_ascending_node_degrees=-124.0)


.. autoclass:: skplatform.satellite.SatelliteMolniya
    :members:
    :special-members: __init__


SatelliteSimpleGeostationary
============================
Create a satellite which follows a Geostationary orbit. This model is a very simple, idealized geostationary orbit that stays
fixed over a given geographic point.  It is useful for modelling/simulation work. For example, the following code will
create a geostationary orbit directly above -80 degree lonngitude::

    geostat = SatelliteSimpleGeostationary( -80.0 )

.. autoclass:: skplatform.satellite.SatelliteSimpleGeostationary
    :members:
    :special-members: __init__


.. x SatelliteTudat
.. x --------------
.. x SatelliteTudat is an experimental class and is not yet available for general usage.
.. x
.. x .. autoclass:: skplatform.satellite.satellitetudat.SatelliteTudat
.. x   :members:
.. x    :special-members: __init__

..  _satellitebase_class:

SatelliteBase
=============
The SatelliteBase is the base class that all satellites inherit from. Most users do not use this class.

.. autoclass:: skplatform.satellite.SatelliteBase

    **Platform Locator methods**
    These are methods required to support :class:`~.PlatformLocator`

    .. automethod:: update_position

    **Satellite Abstract methods**

    .. automethod:: update_eci_position

    .. automethod:: period

   **Regular methods**

    .. automethod:: eciposition

    .. automethod:: ecivelocity

    .. automethod:: equator_crossing

    .. automethod:: set_orbit_number_from_last_equator_crossing

    .. automethod:: orbit_number

    .. automethod:: _set_current_state

Sun, Moon, Stars
=================
skplatform provides a small wrapper around astropy to support stars

    ..  autofunction:: skplatform.celestialbodies.star_unitvector_itrf

    .. autofunction:: skplatform.celestialbodies.solsys_body_vector_itrf





