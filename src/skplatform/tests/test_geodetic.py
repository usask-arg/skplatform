import unittest
import numpy as np
from skplatform.geodesy import Geodetic


class GeodeticTests(unittest.TestCase):

    def test_llh_to_xyz_to_llh(self):
        """
        Tests the coordinate conversion back and forth between llh and xyz. Tests ~16 million settings distributed over
        the sphere and over a range of heights. We include both poles which are known special cases. Everything must
        agree to within 7 decimal places. Many agree to much better.
        """
        geo = Geodetic()
        lats = np.linspace(90, -90, 200)
        lngs = np.linspace(0, 360, 200)
        alts = np.linspace(-4000.0, 80000.0, 400)
        latlonalt = np.zeros([lats.size, lngs.size, alts.size, 3])              # Create a large array to hold all the permutations
        for ilat in range(lats.size):
            latlonalt[ilat, :, :, 0] = lats[ilat]                               # set latitude.
        for ilng in range(lngs.size):
            latlonalt[:, ilng, :, 1] = lngs[ilng]                               # set longitude.
        for ihgt in range(alts.size):
            latlonalt[:, :, ihgt, 2] = alts[ihgt]                               # set height.

        xyz = geo.xyz_from_llh(latlonalt)                                       # Convert geodetic to geocentric
        llh = geo.llh_from_xyz(xyz)                                             # and then convert it back.
        diff = llh - latlonalt                                                  # we should get the same answer back
        maxdlat = np.max(np.abs(diff[..., 0:1]))
        maxdlng = np.max(np.abs(diff[..., 1:2]))
        maxdhgt = np.max(np.abs(diff[..., 2:3]))
        self.assertAlmostEqual(maxdlat, 0.0, 7)
        self.assertAlmostEqual(maxdlng, 0.0, 7)
        self.assertAlmostEqual(maxdhgt, 0.0, 7)
#        print("Max lat diff = {}".format(maxdlat))
#        print("Max lng diff = {}".format(maxdlng))
#        print("Max hgt diff = {}".format(maxdhgt))

    # -----------------------------------------------------------------------------
    #               test_tangent_point_techniques
    # -----------------------------------------------------------------------------
    def test_unit_vectors(self):

        geo = Geodetic()
        cos30 = np.cos(np.radians(30.0))
        sin30 = np.sin(np.radians(30.0))
        cos45 = np.cos(np.radians(45.0))
        sin45 = np.sin(np.radians(45.0))
        tests = [((0.000, +0.0000, 500000.0), ((0, 0, 1), (0, 1, 0), (-1, 0, 0))),
                 ((0.000, +90.000, 500000.0), ((0, 0, 1), (-1, 0, 0), (0, -1, 0))),
                 ((0.000, +180.00, 500000.0), ((0, 0, 1), (0, -1, 0), (1, 0, 0))),
                 ((0.000, +270.00, 500000.0), ((0, 0, 1), (1, 0, 0), (0, 1, 0))),
                 ((+30.0, +0.0000, 500000.0), ((-sin30, 0, cos30), (0, 1, 0), (-cos30, 0, -sin30))),
                 ((-30.0, +0.0000, 500000.0), ((+sin30, 0, cos30), (0, 1, 0), (-cos30, 0, +sin30))),
                 ((+30.0, +45.000, 500000.0), ((-sin30 * cos45, -sin30 * sin45, cos30), (-sin45, cos45, 0), (-cos30 * cos45, -sin45 * cos30, -sin30))),
                 ((-30.0, -45.000, 500000.0), ((+sin30 * cos45, -sin30 * sin45, cos30), (+sin45, cos45, 0), (-cos30 * cos45, +sin45 * cos30, +sin30))),
                 ((+90.0, +68.345, 500000.0), ((-1, 0, 0), (0, 1, 0), (0, 0, -1))),
                 ((-90.0, 0.00000, 500000.0), ((+1, 0, 0), (0, 1, 0), (0, 0, 1))),
                 ]
        for entry in tests:
            location = np.array(entry[0])
            answer = np.array(entry[1])
            north, east, down = geo.xyz_north_east_down(llh=location)
            diff1 = np.max(np.abs(north - answer[0, :]))
            diff2 = np.max(np.abs(east - answer[1, :]))
            diff3 = np.max(np.abs(down - answer[2, :]))
            if np.abs(location[0] != 90.0):
                self.assertAlmostEqual(diff1, 0.0, 7)
                self.assertAlmostEqual(diff2, 0.0, 7)
            self.assertAlmostEqual(diff3, 0.0, 7)

    # -----------------------------------------------------------------------------
    #               test_tangent_point_techniques
    # -----------------------------------------------------------------------------
    def test_tangent_point_techniques(self):
        '''
        Compares the tangent point calculated by the two techniques implement in the class. The two algorithms
        should locate the tangent point within 0.3 m of each other for all geometries tested. We place the observer
        at 500 km at 0 degrees longitude ranging from a latitude of -90 to +90 in 25 steps.

        At every observer position we look at 1000 elevations and 400 azimuths. The elevation ranges from 20.7 to 21.7
        degrees (which gives a height range from around 11 km to 58 km) and the look vector azimuth is uniformly spread
        over the 360 degrees of compass bearing. All told, we compare around 10 million tangent point calculations.

        :return:
        '''

        geo = Geodetic()
        obsvlat = np.linspace(-90.0, 90.0, 25)
        elevate = np.linspace(20.7, 21.7, 1000)
        bearing = np.linspace(0.0, 360.0, 400)
        Nlat = obsvlat.size
        Nelev = elevate.size
        Nazi = bearing.size
        obsllh = np.empty([Nlat, Nelev, Nazi, 3])
        for i in range(Nlat):
            obsllh[i, :, :, 0] = obsvlat[i]
        obsllh[:, :, :, 1] = 0.0
        obsllh[:, :, :, 2] = 500000.0

        obsxyz = geo.xyz_from_llh(obsllh)
        north, east, down = geo.xyz_north_east_down(llh=obsllh)
        cosbear = np.cos(np.radians(bearing)[..., np.newaxis])              # Add the extra axis so we can broadcast through last dim
        sinbear = np.sin(np.radians(bearing)[..., np.newaxis])
        coselev = np.cos(np.radians(elevate)[..., np.newaxis, np.newaxis])  # add the extra 2 axes so we can broadcast through the last 2 dimensions
        sinelev = np.sin(np.radians(elevate)[..., np.newaxis, np.newaxis])

        horiz = north * cosbear + east * sinbear                            # Calculate the horizontal component of the look vector, broadcasting the sin and cos of the bearing
        lookxyz = horiz * coselev + down * sinelev                          # Calculate the look vector, broadcasting the sin and cos of the elevation

        tanpt_geocentric1 = geo.xyz_tangent_point_location(obsxyz, lookxyz)
        tanpt_geocentric2 = geo.xyz_tangent_point_location_legacy(obsxyz, lookxyz)
        diff = np.abs(tanpt_geocentric1 - tanpt_geocentric2)
        delta = np.linalg.norm(diff, axis=-1)
        maxdelta = np.max(delta)
        self.assertLess(maxdelta, 0.30)
#        tanpt_llh = geo.llh_from_xyz(tanpt_geocentric2)
#        mintanheight = np.min(tanpt_llh[..., 2:3])
#        maxtanheight = np.max(tanpt_llh[..., 2:3])

    # -----------------------------------------------------------------------------
    #               test_shell_intercepts
    # -----------------------------------------------------------------------------
    def test_shell_intercepts(self):
        '''
        Compares the altitude intercepts calculated by this version of Geodetic with values calculated with the legacy C++
        version distributed with sasktran. It calculates the intersection location for 20 lines of sight and compares it to the
        same calculation using the legacy code. They should generally agree to better than 0.05 meters and we test for
        0.1 meters.
        '''

        geo = Geodetic()
        Nlook = 20
        H = np.linspace(15000.0, 80000.0, Nlook)
        lookangle = np.zeros([Nlook]) + 21.7
        obsllh = np.zeros([Nlook, 3])
        obsllh[..., 0:1] = np.linspace(-89, 89, Nlook)[:, np.newaxis]
        obsllh[..., 1:2] = 0.0
        obsllh[..., 2:3] = 500000.0

        obsxyz = geo.xyz_from_llh(obsllh)
        north, east, down = geo.xyz_north_east_down(llh=obsllh)
        lookxyz = north * np.cos(np.radians(lookangle)[..., np.newaxis]) + down * np.sin(np.radians(lookangle)[..., np.newaxis])
        shellentrance, shellexit, status = geo.xyz_altitude_intercepts(H, obsxyz, lookxyz)

        skshellentrance = np.array([
            [2268471.02325792, 0.00000000, -5957149.92325174],
            [3148326.68127529, 0.00000000, -5549507.94742907],
            [3971150.74737441, 0.00000000, -5002163.51189557],
            [4706059.59884038, 0.00000000, -4327909.59135542],
            [5329806.46891977, 0.00000000, -3543484.05004515],
            [5824136.84041809, 0.00000000, -2668871.16479160],
            [6175314.25308364, 0.00000000, -1726440.05355786],
            [6373975.23263993, 0.00000000, -740096.29872236],
            [6414977.66786679, 0.00000000, 265479.74932607],
            [6297197.96036337, 0.00000000, 1265469.81547996],
            [6023307.62107467, 0.00000000, 2235428.99006860],
            [5599567.76928594, 0.00000000, 3151702.96974368],
            [5035660.56208531, 0.00000000, 3991800.31262670],
            [4344549.71714516, 0.00000000, 4734768.91897148],
            [3542338.87958080, 0.00000000, 5361612.42908031],
            [2648085.07684005, 0.00000000, 5855758.80149096],
            [1683520.26974417, 0.00000000, 6203570.85839508],
            [672645.56330515, 0.00000000, 6394864.18781071],
            [-358821.18177218, 0.00000000, 6423383.08825766],
            [-1384140.02031840, 0.00000000, 6287177.73544280],
        ])
        skshellexit = np.array([
            [2676942.27428444, 0.00000000, -5786282.69761792],
            [3652888.98540664, 0.00000000, -5233384.00026220],
            [4502645.60660385, 0.00000000, -4532980.02472563],
            [5212899.02831040, 0.00000000, -3706346.48748189],
            [5769085.49220013, 0.00000000, -2777275.87874774],
            [6158868.48190297, 0.00000000, -1771581.09939123],
            [6373221.42741500, 0.00000000, -716716.86878932],
            [6406988.47098928, 0.00000000, 358620.99847029],
            [6259291.83990892, 0.00000000, 1424977.69150543],
            [5933853.76850256, 0.00000000, 2452744.80583541],
            [5439208.75092746, 0.00000000, 3412903.17243315],
            [4788754.87304300, 0.00000000, 4277884.04099939],
            [4000596.09079895, 0.00000000, 5022500.29846643],
            [3097145.90771167, 0.00000000, 5624879.22558630],
            [2104488.81527875, 0.00000000, 6067318.61080514],
            [1051525.18832378, 0.00000000, 6336988.62355512],
            [-31048.77355950, 0.00000000, 6426414.65503854],
            [-1111849.21302771, 0.00000000, 6333697.08947399],
            [-2159914.95066628, 0.00000000, 6062456.60543426],
            [-3145734.50571867, 0.00000000, 5621526.64889211],
        ])

        diff1 = np.max(np.abs(skshellentrance - shellentrance))
        diff2 = np.max(np.abs(skshellexit - shellexit))
        self.assertLess(diff1, 0.1)
        self.assertLess(diff2, 0.1)

        # The above variables were generated by running the code below
        # if skshellentrance is None:
        #     import sasktran as sk
        #     skgeo = sk.Geodetic()
        #     skshellentrance = np.empty_like(shellentrance)
        #     skshellexit = np.empty_like(shellexit)
        #     for i in range(Nlook):
        #         entrancept, exitpt = skgeo.altitude_intercepts(H[i], obsxyz[i, :], lookxyz[i, :])
        #         skshellentrance[i,:] = entrancept
        #         skshellexit[i, :]= exitpt
        #
        #     print('skshellentrance = np.array([')
        #     for i in range(Nlook):
        #         print('[{:18.8f}, {:18.8f}, {:18.8f}],'.format(skshellentrance[i, 0], skshellentrance[i, 1], skshellentrance[i, 2]))
        #     print('])')
        #
        #     print('skshellexit = np.array([')
        #     for i in range(Nlook):
        #         print('[{:18.8f}, {:18.8f}, {:18.8f}],'.format(skshellexit[i, 0], skshellexit[i, 1], skshellexit[i, 2]))
        #     print('])')

    def test_look_at_tangent_point(self):
        '''
        Test the code that finds a look vector that looks at a specific tangent altitude.
        '''

        geo = Geodetic()
        obsllh = np.array([67.0, 0.0, 500000.0])                                                # Position the satellite at a given location, lat,lon, height
        obsxyz = geo.xyz_from_llh(obsllh)                                                       # and get the geocentric xyz value
        tanH = 30000.0                                                                          # Set the target tangent altitude
        north, east, down = geo.xyz_north_east_down(llh=obsllh)                                 # Get the local geographic unit vectors
        azi = np.radians(56.0)                                                                  # set the azimuth o fthe look vector
        lookxyz = north * np.cos(azi) + east * np.sin(azi)                                      # and setup the horizontal look vector direction
        lookv = geo.xyz_lookvector_from_tangent_altitude(tanH, obsxyz, lookxyz)                 # then get the look vector in this plane towards the target tangent altitude

        # sklook = skgeo.from_tangent_altitude(tanH, obsxyz, lookxyz)                           # the following numbers were generated by running the following line with
        sklook = np.array([-0.6219514733030740, 0.7721580709962425, -0.1301855531597018])       # the sasktran C++ version

        tanpt_geocentric = geo.xyz_tangent_point_location(obsxyz, lookv)                        # As a check, Locate the tangent point using the lookv
        tanptllh = geo.llh_from_xyz(tanpt_geocentric)                                           # and get the lat/lon height

        diff1 = np.linalg.norm(sklook - lookv)                                                  # get the difference between the python geodetic and the C++/sasktran geodetic calculation
        diff2 = np.abs(tanptllh[2] - tanH)                                                      # and get the altitude offset of the target tangent altitude form the calculated atitude.
        self.assertAlmostEqual(diff1, 0.0, 5)                                                   # assert both answers are almost equal. It seems that 6 decimal places is not feasible
        self.assertLess(diff2, 0.1)                                                             # and make sure we are close to our target tangent altitude


if __name__ == "__main__":
    tests = GeodeticTests()
    unittest.main()
