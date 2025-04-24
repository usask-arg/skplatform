####################
skplatform.scripting
####################

The ``scripting`` sub-module is a framework built on top of the main package of
``skplatform`` with the goal of providing a simple interface for simulating the
orbit of a satellite and the measurement locations of an on-board instrument.
It provides a scripting interface, a command line interface and classes for
use in the ``hawcsimulator``.

Usage
=====

API and Command Line Interface
------------------------------

``skplatform.scripting`` provides a scripting api and a command line interface.
The command line can be used simply::

    python cli.py --instrument aos-sky --end "2019-11-01 00:00:01" --output "~/data/aos/" --threads 16

The help for the CLI can be accessed using::

    python cli.py --help

which will print out a list of all the arguments with detailed descriptions.
The CLI has two main limitations. It does not support custom instruments and
it does not support non-uniform frequencies.

The CLI is simply a wrapper around the scripting interface, the function:
``skplatform.scripting.interface.simulate_orbit``. ``simulate_orbit()``
supports custom instruments in addition to the pre-defined ones.
Additionally, non-uniform frequencies can be used. If your IDE does not
automatically provide help for the function, using::

    from skplatform.scripting.interface import simulate_orbit
    help(simulate_orbit)

will print out a list of inputs with detailed descriptions.

Refer to ``src/skplatform/scripting/examples`` for more examples.

Concepts
========

Observations and Measurements
-----------------------------

During a simulation the satellite position is calculated and stored based on
the instrument frequency. In the context of this submodule, each time the
satellite position is calculated is termed an *observation*. At each
observation the lines-of-sight for a theoretical instrument are calculated.
For a limb instrument the tangent point of each LOS is determined and saved.
For a nadir instrument the location of the ground intercept of each LOS is
determined and saved. For brevity, the tangent point/ground intercept of each
LOS is called a *measurement*.

Observation Frequency
---------------------

The temporal spacing of observations is termed *frequency*. In the simple
case when frequency is a single number, it is the frequency at which
observations are made, ie. the number of seconds between observations.
However, the term is misused somewhat, in that the frequency can be an array
specifying a set of non-uniform intervals. The pattern of observations
described in the frequency array is repeated over the length of the
simulation. Each entry in the array, except the last, sets the time interval,
in seconds, between observations. So the first observation occurs immediately,
the next observation occurs ``frequency[0]`` seconds after the first, and so
forth. Generally, the interval between observations *n-1* and *n* is
``frequency[n-1]``. The last entry in the frequency array is the time interval
before the set is repeated, ie. it is the interval between the last
observation of the current set and the first observation of the next set.
Whether the frequency is a single number or an array, the simulation can
end between observations.

Instrument Model
----------------

The instrument in ``skplatform.scripting`` is abstract, and in both limb and
nadir orientations it is modelled as a rectangular grid of theoretical pixels.
Each pixel is associated with a single line-of-sight (LOS). The center of the
grid is the boresight, and is considered the direction of the instrument. All
other LOS are defined as being at an angle relative to the boresight. Note
that the boresight is also an abstract concept, and instruments with an even
number of pixels on either or both axes will not have LOS calculated for the
boresight. For nadir instruments the boresight cannot be configured and always
looks straight down from the satellite. For limb instruments the boresight has
a target altitude and angle.

Instrument Classes
------------------

Support for experimentation with different instrument configurations is
provided by the ``LimbInstrument`` and ``NadirInstrument`` classes. Simply,
import the desired class, instantiate it::

    from skplatform.scripting.instruments import LimbInstrument
    inst = LimbInstrument(tle_lines,
                        along_fov=0.50, across_fov=7.0,
                        along_pixels=5, across_pixels=30,
                        frequency=np.array([2.0, 5.0, 9.0]),
                        look_angle_deg=180.0,
                        target_altitude_m=20_000.0)

and pass it to ``simulate_orbit()``::

    ds = simulate_orbit(instrument=inst, output=None, tle=tle, start=start, end=end, threads=1)


If more flexibility is needed custom Instrument classes can also be derived
from the ``skplatform.scripting.instrument.PropagatorInstrument`` class.
Refer to the source code for the details.
