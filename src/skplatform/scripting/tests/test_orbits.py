

def default_tle():
    """
    Returns a TLE set as a list of two strings.
    """
    # taken from AOS_Sky_Descending_v2025-01-29_Epoch2019.tle
    # directly entered so don't rely on the file
    return ['1 99999U 99999 A  19213.00000000  .00000000  00000-0  00000-0 0    10\n',
            '2 99999 097.2222 331.8182 0001962 091.0782 359.0176 15.37843164    16\n']


def setup_sat():
    """
    Returns a Skplatform SatelliteSGP4 instance. The satellite is initialized with the default TLE.
    """
    from skplatform.satellite import SatelliteSGP4
    sat = SatelliteSGP4(twolines=default_tle())
    sat.use_polar_motion = True
    return sat


def setup_limb_propagator():
    """
    Returns a LimbPropagator instance, with observations every 100s for 7.5 hours (approximately 5 orbits).
    """
    import numpy as np
    from datetime import timedelta
    from skplatform import Platform
    from skplatform.scripting.propagator import LimbPropagator
    from skplatform.scripting.rotations import calculate_2d_rotation_matrices

    sat = setup_sat()
    platform = Platform(platform_locator=sat)
    # approximately 5 orbits, with measurements every 100s
    times = np.array([timedelta(seconds=dt) for dt in range(0, 90*60*5, 100)])
    rot_mat = calculate_2d_rotation_matrices(0, 5.0, 1, 64)
    boresight = 180.0
    altitude = 20_000.0
    limb_propagator = LimbPropagator(platform, times, rot_mat, boresight, altitude)
    return limb_propagator


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
    import numpy as np
    return np.nanmax(np.abs( (a-b)))


def test_eci():
    """
    Compares the ECI position of the satellite calculated using skplatform, Astropy and Skyfield.
    Comparisons are between newly calculated skplatform positions and stored results in orbit_test_data.py.
    """
    import numpy as np
    from numpy.testing import assert_allclose, assert_almost_equal

    from skplatform import Platform
    import skplatform.scripting.tests.orbit_test_data as orbit_test_data

    scn = setup_limb_propagator()

    # get the skplatform velocity and positions
    sat = setup_sat()
    sk_plat = Platform(platform_locator=sat)

    # store the eci coordinates
    sat_eci_pos = np.zeros((scn.num_measurements, 3))

    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]
        # SkPlatform doesn't store the eci values
        sk_plat.platform_locator.update_eci_position(time)

        sat_eci_pos[t, :] = sk_plat.platform_locator.eciposition()

    # From orbit_comparisons.py
    # Comparison of ECI (TEME) Coordinates
    #     skplatform vs Skyfield:         4.34229e-08 m   2.56366e-10%
    #     skplatform vs sk->Astropy:      0.00000e+00 m   0.00000e+00%
    #     skplatform vs Skyfield Sat:     9.96515e-08 m   1.64145e-10%
    #     skplatform vs SGP4 Sat:         7.63685e-08 m   1.64145e-10%
    # # since we are saving to disk and reading back values a small amount of precision is lost

    # vs_self = (relative_difference(sat_eci_pos.flatten(), orbit_arrays.sat_eci.flatten()))
    sky_rel = (absolute_difference(sat_eci_pos.flatten(), orbit_test_data.sky_eci.flatten()))
    astro_rel = (absolute_difference(sat_eci_pos.flatten(), orbit_test_data.astro_eci.flatten()))
    sky_sat_rel = (absolute_difference(sat_eci_pos.flatten(), orbit_test_data.sky_sat_eci.flatten()))
    astro_sat_rel = (absolute_difference(sat_eci_pos.flatten(), orbit_test_data.astro_sat_eci.flatten()))

    assert_allclose(sat_eci_pos.flatten(), orbit_test_data.sat_eci.flatten(), err_msg='Skplatform ECI does not match itself')
    assert_almost_equal(sky_rel, 4.34229e-08, err_msg='Skplatform ECI does not match ECI converted using Skyfield')
    assert_almost_equal(astro_rel, 0.00000e+00, err_msg='Skplatform ECI does not match ECI converted from Astropy TEME ')
    assert_almost_equal(sky_sat_rel, 9.96515e-08, err_msg='Skplatform ECI does not match Skyfield Satellite ECI calculated from time')
    assert_almost_equal(astro_sat_rel, 7.63685e-08, err_msg='Skplatform ECI does not match Astropy Satellite ECI calculated from time')


def test_ecef():
    """
    Compares the ECEF (ITRS) position of the satellite calculated using skplatform, Astropy and Skyfield.
    Comparisons are between newly calculated skplatform positions and stored results in orbit_test_data.py.
    """
    import numpy as np
    from numpy.testing import assert_allclose, assert_almost_equal

    from skplatform import Platform
    import skplatform.scripting.tests.orbit_test_data as orbit_test_data

    scn = setup_limb_propagator()

    sat = setup_sat()
    sk_plat = Platform(platform_locator=sat)

    sat_ecef = np.zeros((scn.num_measurements, 3))
    for t in range(scn.num_measurements):
        time = scn.measurement_times[t]
        # SkPlatform doesn't store the eci values
        sk_plat.platform_locator.update_eci_position(time)
        sat_ecef[t, :] = sk_plat.platform_locator.position

    # ================================================================================
    # check_eci_to_itrs()
    # Checking ECI to ITRS conversion
    # ================================================================================
    # Comparison of ITRS Coordinates
    #         skplatform vs Skyfield: 4.90297e-03 m   8.55150e-05%
    #         skplatform vs Astropy:  1.29157e-05 m   2.22371e-07%
    sky_rel = absolute_difference(sat_ecef.flatten(), orbit_test_data.sky_ecef.flatten())
    astro_rel = absolute_difference(sat_ecef.flatten(), orbit_test_data.astro_ecef.flatten())
    assert_allclose(sat_ecef.flatten(), orbit_test_data.sat_ecef.flatten(), err_msg='Skplatform ECEF does not match '
                                                                                    'itself')
    assert_almost_equal(astro_rel, 1.29157e-05, err_msg='Skplatform ECEF does not match Astropy calculated ECEF')
    assert_almost_equal(sky_rel, 4.90297e-03, err_msg='Skplatform ECEF does not match Skyfield calculated ECEF')


def test_itrs_to_lla():
    """
    Compares the latitude, longitude and altitude of the satellite calculated using skplatform, Astropy and Skyfield.
    The input is the ECEF (ITRS) position calculated by skplatform.
    Comparisons are between newly calculated skplatform values and stored results in orbit_test_data.py.
    """
    import numpy as np
    from numpy.testing import assert_allclose, assert_almost_equal

    import skplatform.scripting.tests.orbit_test_data as orbit_test_data

    scn = setup_limb_propagator()

    sk_lla = np.zeros((scn.num_measurements, 3))

    scn.propagate()

    sk_lla[:, 0] = scn.ds.observer_latitude.values
    sk_lla[:, 1] = scn.ds.observer_longitude.values
    sk_lla[:, 2] = scn.ds.observer_altitude.values

    # From orbit_comparisons.py
    # skplatform vs Astropy:
    #     Lat: 0.00000e+00 m      0.00000e+00%
    #     Lon: 5.68434e-14 m      8.22143e-13%
    #     Alt: 0.00000e+00 m      0.00000e+00%
    # skplatform vs Skyfield:
    #     Lat: 1.94289e-08 m      1.52655e-07%
    #     Lon: 2.50111e-12 m      2.68929e-10%
    #     Alt: 4.16185e-08 m      9.09391e-12%
    astro_rel = (absolute_difference(sk_lla.flatten(), orbit_test_data.astro_lla.flatten()))
    sky_rel = (absolute_difference(sk_lla.flatten(), orbit_test_data.sky_lla.flatten()))
    assert_allclose(sk_lla.flatten(), orbit_test_data.sk_lla.flatten(), err_msg='Skplatform LLA from ECEF does not match itself')
    assert_almost_equal(astro_rel.flatten(), 5.68434e-14, err_msg='Skplatform LLA from ECEF does not match LLA converted from ECEF using Astropy')
    assert_almost_equal(sky_rel.flatten(), 4.16185e-08, err_msg='Skplatform LLA from ECEF does not match LLA converted from ECEF using Skyfield')


def test_orbit_to_lla():
    """
    Compares the latitude, longitude and altitude of the satellite calculated using skplatform, Astropy and Skyfield.
    The values are calculated from the observation times using each library's propagator interface. The results are
    the other libraries are independent of skplatform.
    Comparisons are between newly calculated skplatform values and stored results in orbit_test_data.py.
    """
    from numpy.testing import assert_allclose, assert_almost_equal

    import skplatform.scripting.tests.orbit_test_data as orbit_test_data

    scn = setup_limb_propagator()
    scn.propagate()

    assert_allclose(scn.ds.observer_latitude.values.flatten(), orbit_test_data.sk_lat, err_msg='Skplatform latitude does not match itself')
    assert_allclose(scn.ds.observer_longitude.values.flatten(), orbit_test_data.sk_lon, err_msg='Skplatform longitude does not match itself')
    assert_allclose(scn.ds.observer_altitude.values.flatten(), orbit_test_data.sk_alt, err_msg='Skplatform altitude does not match itself')

    # From orbit_comparisons.py
    # Skplatform vs SGP4->Astropy
    #         Lat: 8.02913e-13 deg    1.64193e-10%
    #         Lon: 1.10674e-10 deg    2.22317e-07%
    #         Alt: 5.35510e-09 m      1.17375e-12%
    astro_rel_lat = absolute_difference(scn.ds.observer_latitude.values.flatten(), orbit_test_data.astro_lat.flatten())
    astro_rel_lon = absolute_difference(scn.ds.observer_longitude.values.flatten(), orbit_test_data.astro_lon.flatten())
    astro_rel_alt = absolute_difference(scn.ds.observer_altitude.values.flatten(), orbit_test_data.astro_alt.flatten())
    assert_almost_equal(astro_rel_lat, 8.02913e-13, err_msg='Skplatform latitude does not match Astropy calculated latitude')
    assert_almost_equal(astro_rel_lon, 1.10674e-10, err_msg='Skplatform longitude does not match Astropy calculated longitude')
    assert_almost_equal(astro_rel_alt, 5.35510e-09, err_msg='Skplatform altitude does not match Astropy calculated altitude')

    # From orbit_comparisons.py
    # Skplatform vs Skyfield
    #         Lat: 3.20314e-08 deg    1.19172e-05%
    #         Lon: 1.39857e-07 deg    8.55053e-05%
    #         Alt: 4.76703e-06 m      1.02396e-09%
    sky_rel_lat = absolute_difference(scn.ds.observer_latitude.values.flatten(), orbit_test_data.sky_lat.flatten())
    sky_rel_lon = absolute_difference(scn.ds.observer_longitude.values.flatten(), orbit_test_data.sky_lon.flatten())
    sky_rel_alt = absolute_difference(scn.ds.observer_altitude.values.flatten(), orbit_test_data.sky_alt.flatten())
    assert_almost_equal(sky_rel_lat, 3.20314e-08, err_msg='Skplatform latitude does not match Skyfield calculated latitude')
    assert_almost_equal(sky_rel_lon, 1.39857e-07, err_msg='Skplatform longitude does not match Skyfield calculated longitude')
    assert_almost_equal(sky_rel_alt, 4.76703e-06, err_msg='Skplatform altitude does not match Skyfield calculated altitude')


def test_planar():
    """
    Test how planar the orbit is. An ideal orbit would reside entirely on a single plane when plotted in the Earth
    Centered Inertial (ECI) frame. The SGP4 propagator includes perturbations that prevent the orbit from being
    perfectly planar, but a significant deviance is not expected. Reference values were not found, but values of
    0.3 degrees deviation were calculated. A perfectly planar orbit or an orbit with significant deviance indicate
    a problem.

    Planarity is checked by comparing the angle between the cross-product of position and velocity of each observation
    to that of every other observation.
    """
    import numpy as np
    from numpy.testing import assert_allclose, assert_almost_equal

    from skplatform import Platform
    import skplatform.scripting.tests.orbit_test_data as orbit_test_data

    scn = setup_limb_propagator()

    sat = setup_sat()
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

    # From orbit_comparisons.py
    # Difference 1.38805e-03  3.01884e-01 deg Mean: 2.37706e-04       Median: 1.28888e-04
    assert_allclose(np.array(compare).flatten(), orbit_test_data.planar.flatten(), atol=1e-15,
                    err_msg='The deviations from a true planar orbit do not match previous results')


if __name__ == "__main__":
    test_eci()
    print('Tests comparing ECI passed')

    test_ecef()
    print('Tests comparing ECEF passed')

    test_itrs_to_lla()
    print('Tests comparing LLA from ECEF passed')

    test_orbit_to_lla()
    print('Tests comparing LLA from propagators passed')

    test_planar()
    print('Tests checking orbit planarity passed')
