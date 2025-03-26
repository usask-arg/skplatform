..  _highlevelif:

###############################
High-Level Simulation Interface
###############################
The ``skplatform`` package provides support to simulate platforms such as spacecraft, aircraft, balloons and
ground-based sites. Instruments are mounted on platforms at given orientations, similar to how physical instruments are
bolted on to real platforms, and the platform is positioned and rotated to collect the desired set of measurements.

Its goal within the ``skretrieval`` framework is to assist in the generation of the set of times, positions and look vectors.
Users who have their own techniques for calculating positions and look vectors for the forward model can continue using
their own code or, if they wish, switch over to using the :class:`~.Platform` class.

Class :class:`~.Platform` provides access to two primary aspects,

#. Simulation of various platforms used in remote atmospheric sensing.
#. Generation of sets of platform position and orientation using various positioning and pointing techniques.

The :class:`~.Platform` class utilizes several :ref:`controlframes` and provides rotation matrices to transform between the
frames. The user will typically define instrument lines of sight in the :ref:`icf`. The instrument is optionally
mounted on a platform and the platform rotated in space and time to execute the desired measurement plan. Users will normally
transform the instrument look vectors to the :ref:`ecef` for input to radiative transfer codes but the platform can also
transform the look vector to the local geodetic control frame or platform control frame if preferred.

Class :class:`~.Platform` provides a high level interface which can be used to calculate platform positions and
orientations needed to simulate a set of measurements. The class, as provided out-of-the-box, provides access to a library
of common techniques that that should meet most user needs but the API allows new techniques to be developed and added
by users. The problem usually reduces into two parts,

#. selection of an instrument positioning technique.
#. selection of an instrument pointing technique.

Example
=======
Consider a simulation where an instrument onboard a sun-synchronous satellite at an altitude of 600 km needs to look in the
limb at a tangent point at 35 km altitude over Saskatoon, Saskatchewan, located at (52.0N, -107.0E). The observation
will be made at approximately 2020-07-15 15:00 UTC. The satellite is configured so the instrument is always looking
at a bearing of 45 degrees from local North at the satellite location.

We shall demonstrate three techniques that could be used to *solve* the problem.  The three techniques are based
upon significantly different algorithms but, from the user's perspective, are implemented in an almost identical manner.

Option 1: llh
-------------
Use the :class:`~.Platform` class with no special configuration. Set the time and position of the
satellite with the :ref:`llh` positioning technique.  Then configure the platform orientation with the
:ref:`tangent_altitude` pointing technique using :ref:`rollcontrol_limb` roll control so it is looking at that tangent
point. This will generate the required position and look vector for the given conditions. The option has the limitation
that it is difficult to choose the position of the platform so the tangent point is over Saskatoon::

    from skplatform import Platform

    platform = Platform()
    utc = ['2020-07-15T15:00:00.0000000']
    observer = [52.0, -107.0, 600000.0]
    tanpoint = [35000.0, 45.0]
    platform.add_measurement_set(utc, ('llh', observer), ('tangent_altitude', 'limb', tanpoint))
    opticalmeasurements = platform.make_optical_geometry()

If you want to print out the location of the satellite and the tangent point you can use the following snippet of code::

    import sasktran as sk

    geo = sk.Geodetic()
    entry = opticalmeasurements[0]
    geo.from_xyz(entry.observer)
    obslat = geo.latitude
    obslng = geo.longitude
    obshgt = geo.altitude / 1000.0
    geo.from_tangent_point(entry.observer, entry.look_vector)
    tanlat = geo.latitude
    tanlng = geo.longitude
    tanhgt = geo.altitude / 1000.0
    print('Satellite location = ({:5.2f}N,{:6.2f}E) at a height of {:6.2f} km'.format(obslat, obslng, obshgt))
    print('Tangent location   = ({:5.2f}N,{:6.2f}E) at a height of {:6.2f} km'.format(tanlat, tanlng, tanhgt))

and you should get::

    Satellite location = (52.00N,253.00E) at a height of 600.00 km
    Tangent location   = (63.65N,291.78E) at a height of  35.00 km

The tangent point is at 35 km as requested but it is not over Saskatoon, which is not surprising as we, perhaps erroneously,
put the satellite directly over Saskatoon. To get the tangent point directly over Saskatoon you could iteratively adjust the satellite
position until you find teh right locations. But there is an easier way demonstrated below in option 2.

Option 2: looking_at_llh
------------------------
The second option provides a more sophisticated technique to position the satellite so the resultant tangent
point will end up above Saskatoon. In this case set the time and position of the satellite are set with the :ref:`looking_at_llh`
positioning technique; this technique has 5 parameters that position the satellite at the requested altitude so it will
have the requested tangent altitude at the requested location with the requested bearing. It then configures the platform
orientation with the :ref:`tangent_altitude` pointing technique using :ref:`rollcontrol_limb` roll control so it is looking at that tangent
point::

    from skplatform import Platform

    platform = Platform()
    utc = ['2020-07-15T15:00:00.0000000']
    observer = [52.0, -107.0, 35000.0, 45.0, 600000.0]
    tanpoint = [35000.0, 45.0]
    platform.add_measurement_set(utc, ('looking_at_llh', observer), ('tangent_altitude', 'limb', tanpoint))
    opticalmeasurements = platform.make_optical_geometry()

Examination of the location of the satellite and tangent point using the same code as in option 1 above gives::

    Satellite location = (38.23N,226.13E) at a height of 600.00 km
    Tangent location = (52.02N,252.98E) at a height of  35.00 km

Now the tangent point is above Saskatoon and the satellite has been placed at a location so it is looking toward Saskatoon
at a tangent altitude of 35 km and a bearing of 45 degrees. Note that the tangent point is not "*exactly*" above Saskatoon.
It is off by approx 0.02 degrees which is due to limitations in the algorithm used.

Option 3: satellite
-------------------
This option configures the :class:`~.Platform` class with a sun-synchronous satellite orbit propagator.  the satellite
is *flown* along the orbit track by providing a sequence of 100 universal times at 1 minute intervals. The position of
the platform is now extracted using the :ref:`from_platform` positioning technique while the look vector is still extracted
with the :ref:`tangent_altitude` pointing technique and :ref:`rollcontrol_limb` roll control.
This option requires a bit more setup but allows the user to *fly* the satellite around the Earth and generate other
positions and look vectors::

    mjd0 = asktime.ut_to_mjd('2020-07-15T15:00:00.000000')                                                         # Define the time of the ascending node
    satellite = SatelliteSunSync(mjd0,                                                                             # Create a sun-synchronous satellite
                                 orbittype='sgp4',
                                 period_from_altitude=600000.0,
                                 localtime_of_ascending_node_hours=2.75)
    platform = Platform(platform_locator=satellite)                                                                 # Create a platform which can use the sun-synchronous satellite for its position

    utc = mjd0 + np.arange(0, 100) / 1440.0                                                                         # Get measurements every minute along the orbit starting at the ascending node.
    looktanalt = [(35000.0, 45)]                                                                                    # Look at a tangent altitude of 35 km, at a geographic bearing of 45 degrees at the satellite location. Note, the parameter set will be expanded to N measurements
    platform.add_measurement_set(utc, ('from_platform',), ('tangent_altitude', 'limb', looktanalt))                 # Add the measurement set
    opticalmeasurements = platform.make_optical_geometry()

This technique has the same problem as option 1 in that it is difficult to pre-configure the satellite to fly exactly over Saskatoon
with the proper geometry.  However, we fiddled with the local time of the ascending node of the satellite and if you look at the
38'th entry of the 100 returned ``opticalmeasurement`` values you should see that it has a tangent point close to Saskatoon::

    38 Satellite location = (38.26N,235.59E) at a height of 605.39 km
    38 Tangent location   = (52.09N,262.62E) at a height of  35.00 km

------------------------------------------------------------------------------------------------------------------------

Measurement Sets
================
The user must define what constitutes a set of measurements. This may be as simple as one exposure or something much more
complicated such as a spatial scan of a specific atmospheric region. In our experience, we almost always break
the definition of measurements sets into the following three steps,

#. Define (or acquire) a set of universal times at which the measurements are made.
#. Specify (or acquire) the position of the platform associated with each measurement at each universal time.
#. Specify (or acquire) the pointing of the platform associated with each measurement at each universal time.

Example
-------
Configure a satellite to acquire a height profile of spectra over Saskatoon (52N, -107E) in the summer time from which
atmospheric aerosols can be inferred over Saskatooon .

**Option 1**:  A simple emulator fixes a satellite ~20 degrees south of Saskatoon at 600 km looking northwards in the limb. There
will be 50 measurements scanning from 0 km to 50 km all made at the same location and same time UTC (2020-06-21T18:00:00.000000)

**Option 2:** A sun-synchronous satellite with an ascending node at 12:00LT will fly over Saskatoon on 2020-06-21. It will collect 50 measurements
scanning from 0 to 50 km. The measurements will be one second apart.

The ``skplatform`` package can help with both of these scenarios and the user would choose which scenario they
wish to model.

Internal Representation
-----------------------
Measurement sets are internally stored within class :class:`~.Platform` using an instance of :class:`PositionAndOrientationArray`. This
class stores a set of times, position and rotation matrices used for the measurements. This information is independent of the
instrument, apart from its mounting orientation, and cannot be used by retrieval and radiative transfer codes until it is
converted into a set of times, positions and instrument look vectors, such as those found in :class:`~.OpticalGeometry`.




