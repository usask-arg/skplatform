import numbers
from typing import Union, List, Tuple
import math
import numpy as np
import scipy.optimize


def _unitvector(x: np.ndarray) -> np.ndarray:
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def _unitvector_debugversion(x: np.ndarray) -> np.ndarray:
    mag = np.linalg.norm(x, axis=-1, keepdims=True)
    bad = (mag < 1.0E-15)
    if np.any(bad):
        mag[bad] = 1.0
    val = x / mag
    return val


def _dotproduct(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.nansum(a * b, axis=-1, keepdims=True)


def _vector_is_zero(x: np.ndarray):
    return np.linalg.norm(x, axis=-1) == 0


# -----------------------------------------------------------------------------
#               _HeightOffsetEvaluator
# -----------------------------------------------------------------------------
class _HeightOffsetEvaluator:
    '''
    An internal class used while locating a shell height location. The zbrent function calls the
    class function operator() to evaluate the offset in height of a given point from the target altitude.
    '''

    def __init__(self, geoid: 'Geodetic', tanpoint: np.ndarray, look: np.ndarray, H: float):

        self.geoid = geoid                  # The geoid used for calculations
        self.tangentpoint = tanpoint        # The location of the tangent point of this look vector
        self.look = look                    # The look direction as a unit vector (from the satellite).
        self.location = tanpoint            # Upon exit this will have the location of the height offset.
        self.H = H

    def evaluate(self, x: float):
        self.location = self.tangentpoint - self.look * x
        llh = self.geoid.llh_from_xyz(self.location)
        return llh[2] - self.H

    def __call__(self, x: float):
        return self.evaluate(x)

    def find_crossing_point(self, lmin: float, lmax: float) -> Tuple[np.ndarray, bool]:

        delta = 0.05 * lmin
        answer = self.evaluate(lmin)  # Get the height offset closest to tangent point
        while (answer >= 0.0):
            lmin -= delta
            answer = self.evaluate(lmin)

        answer = self.evaluate(lmax)
        while (answer <= 0.0):
            lmax += 0.05 * lmax
            answer = self.evaluate(lmax)

        answer, r = scipy.optimize.brentq(self, lmin, lmax, xtol=0.1, maxiter=50, full_output=True)

        if r.converged:
            self.evaluate(answer)
        return self.location, r.converged


# -----------------------------------------------------------------------------
#               Geodetic
# -----------------------------------------------------------------------------
class Geodetic():
    '''
    A class to manage both transform and tangent point calculations using oblate spheroid geometry.
    The transform operations allow conversion between geocentric, xyz, coordinates and geodetic, latitude, longitude
    and height coordinates. We also provide right-handed, geographic-based, unit vectors for any location

    The tangent point calculations include finding the location of a tangent point as well as the entrance and exit
    locations of a ray through altitude shells.

    The class is initialized with the default WGS84 oblate spheroid which should be suitable for most work.

    The class is fully implemented in pure python and has been written allow multiple calculations to be performed using
    array-based arithemetic.
    '''

    def __init__(self):
        self.A: float = 0.0
        self.F: float = 0.0
        self._AEPS2 = 0.0
        self._EC2 = 0.0
        self._EC = 0.0
        self._B = 0.0
        self._E4T = 0.0
        self.set_geoid(geoid='wgs84')

    # -----------------------------------------------------------------------------
    #               set_geoid
    # -----------------------------------------------------------------------------
    def set_geoid(self, A: float=None, F: float=None, geoid: str=None):
        '''
        Defines the shape parameters of the ellipsoid that will be used for future calculations. The shape is defined
        by the semi-major axis, ``A``, in meters and the Earth flattening parameter , ``F``. The most common usage is to
        use a standard, pre-defined ``geoid``. This method only needs to be called if the default oblate spheroid used
        in the constructor, `wgs84`, is not suitable.

        Parameters
        ----------
        geoid: str, optional
            Predefined, standard oblate spheroids:

            * 'geocentric': (6371200.0, 0.0)
            * 'iau1976':    (6378140.0, 0.00335281)
            * 'grs80':      (6378137.0, 1.0/298.257222101)
            * 'merit83':    (6378137.0, 1.0/298.257)
            * 'wgs84':      (6378137.0, 1.0/298.257223563)}

            If it is not set the both ``A`` and ``F`` must be set. ``A`` and ``F`` are ignored in this option is set.

        A: float, optional
            The semi-major axis in meters.  This argument is ignored if ``geoid`` is set. Argument ``F`` must also be set if this argument is set.

        F: float, optional
            The Earth flattening parameter. This argument is ignored if ``geoid`` is set. Argument ``A`` must also be set if this argument is set
        '''
        if (geoid is not None):
            knowngeoids = {'geocentric': (6371200.0, 0.0),
                           'iau1976': (6378140.0, 0.00335281),
                           'grs80': (6378137.0, 1.0 / 298.257222101),
                           'merit83': (6378137.0, 1.0 / 298.257),
                           'wgs84': (6378137.0, 1.0 / 298.257223563)}
            key = geoid.lower()
            g = knowngeoids.get(key)
            if (g is not None):
                A = g[0]
                F = g[1]
            else:
                raise ValueError("Unsupported geoid {}. Try 'wgs84` or 'iau1976'".format(geoid))
        else:
            assert ((A is not None) and (F is not None))

        if (F < 0.0 or F >= 1.0):                                                   # Validate ellipsoid parameters.
            raise ValueError('Invalid Earth flattening parameter')
        if (A <= 0.0):
            raise ValueError('Invalid Earth radius parameter')

        # *  Functions of ellipsoid parameters (with further validation of F).

        self.A = A
        self.F = F
        self._AEPS2 = A * A * 1.0E-32
        self._E2 = (2.0 - F) * F
        self._E4T = self._E2 * self._E2 * 1.5
        self._EC2 = 1.0 - self._E2
        if (self._EC2 <= 0.0):
            raise ValueError('Invalid EC2 value')
        self._EC = math.sqrt(self._EC2)
        self._B = A * self._EC

    # -----------------------------------------------------------------------------
    #               semi_major_axis
    # -----------------------------------------------------------------------------
    @property
    def semi_major_axis(self):
        '''
        Semi-major axis of the current ellipsoid in meters
        '''
        return self.A

    # -----------------------------------------------------------------------------
    #               semi_minor_axis
    # -----------------------------------------------------------------------------
    @property
    def semi_minor_axis(self):
        '''
        Semi-minor axis of the current ellipsoid in meters
        '''
        return self.A * (1.0 - self.F)

    # -----------------------------------------------------------------------------
    #               llh_from_xyz
    # -----------------------------------------------------------------------------
    def llh_from_xyz(self, xyz: np.ndarray) -> np.ndarray:
        '''
        Converts positions specified in xyz geocentric coordinates in meters to the equivalent
        geodetic  (latitude, longitude, height) coordinates. The conversion is non-trivial and this implementation
        uses the algorithm developed by Fukushima. In theory, the algorithm is described as an approximation but in practice it
        is numerically comparable to exact solutions operating, which provide solutions at the nanometer level of precision, but
        struggle to solve the cubic equations with any higher precision.

        Notes
        -----
        The algorithm is based upon the ``GCONV2H`` subroutine provided in file ``xgconv2.txt`` which is provided by Toshio Fukushima, `[1]`_
        on his ResearchGate site.

        The code is implemented to efficiently convert a multi-dimensional array of positions in one pass of the routine.

        ..  [1] Fukushima, T., "Transformation from Cartesian to geodetic coordinates accelerated by Halley's method", J.Geodesy (2006) 79: 689-693

        Parameters
        ----------
        llh: array[..., 3], optional
            An array of position vectors given in geocentric, xyz, meters. The array can be multi-dimensional but the last dimension
            must be used for the (x,y,z) fields and its size must be equal to 3, [0=X, 1=Y, 2=Z].

        Returns
        -------
        array [..., 3]
            Returns the geodetic coordinates of the input geocentric locations. The array is the same size and shape as the input array
            but the last dimensions, which is still size 3, is now used for the geodetic coordinates latitude=0, longitude = 1, height =2.
            latitude is given in degrees, -90.0 to +90.0, longitude is given in degrees, usually 0 to 360 (but this is not rigorously enforced),
            and height is given as meters above the geodetic surface (aka sea level)
        '''

        if type(xyz) is not np.ndarray:
            xyz = np.array(xyz, dtype=float)
        if xyz.dtype.name != 'float64':
            xyz = np.array(xyz, dtype='float')
        X = xyz[..., 0:1]          # Our X Y Z must be along the last axis
        Y = xyz[..., 1:2]
        Z = xyz[..., 2:3]

        P2 = X * X + Y * Y                              # Distance from polar axis squared.
        ELONG = np.arctan2(Y, X)                        # set the longitude
        ABSZ = np.abs(Z)                                # Unsigned z-coordinate.
        arepolar = P2 <= self._AEPS2                    # find out if any points are polar
        P = np.sqrt(P2)                                 # Distance from polar axis.
        S0 = ABSZ / self.A                              # Normalization.
        PN = P / self.A
        ZC = self._EC * S0

        C0 = self._EC * PN                                          # Prepare Newton correction factors.
        C02 = C0 * C0
        C03 = C02 * C0
        S02 = S0 * S0
        S03 = S02 * S0
        A02 = C02 + S02
        A0 = np.sqrt(A02)
        A03 = A02 * A0
        D0 = ZC * A03 + self._E2 * S03
        F0 = PN * A03 - self._E2 * C03

        B0 = self._E4T * S02 * C02 * PN * (A0 - self._EC)           # Prepare Halley correction factor.
        S1 = D0 * F0 - B0 * S0
        CC = self._EC * (F0 * F0 - B0 * C0)

        PHI = np.arctan(S1 / CC)                                    # Evaluate latitude and height.
        S12 = S1 * S1
        CC2 = CC * CC
        HEIGHT = (P * CC + ABSZ * S1 - self.A * np.sqrt(self._EC2 * S12 + CC2)) / np.sqrt(S12 + CC2)

        if np.any(arepolar):                                        # Correct any points that are considered polart
            PHI[arepolar] = math.pi / 2.0
            HEIGHT[arepolar] = ABSZ[arepolar] - self._B

        ELONG = np.degrees(ELONG)
        negative = ELONG < 0.0
        ELONG[negative] += 360.0

        southern = Z < 0.0                                          # find any points which are in the southern hemisphere
        PHI[southern] = -PHI[southern]                              # and restore sign of latitude.
        result = np.empty_like(xyz)
        result[..., 0:1] = np.degrees(PHI)
        result[..., 1:2] = ELONG
        result[..., 2:3] = HEIGHT
        return result

    # -----------------------------------------------------------------------------
    #               xyz_from_llh
    # -----------------------------------------------------------------------------
    def xyz_from_llh(self, llh: Union[np.ndarray, Tuple[float, float, float], List[float]]) -> np.ndarray:
        '''
        Converts positions specified in geodetic coordinates (lat degrees,lon degrees, height meters) to the equivalent
        geocentric (x,y,z meters ) coordinates. The conversion is implemented to it can efficiently process a multi-dimensional
        array of positions in one pass.

        Notes
        -----
        The conversion from geodetic to geocentric is straight-forward and quick. For more details see the reference:
        Page K11 Astronomical Almanac 1990.

        Parameters
        ----------
        llh: array[..., 3], optional
            An array of position vectors given in geodetic latitude, longitude, height (degrees, degrees, meters). The
            array can be multi-dimensional but the last dimension is reserved for the lat, lon and height fields and it size
            must be equal to 3, [0=lat, 1=long, 2=height].

        Returns
        -------
        array[..., 3]:
            The geocentric, xyz, locations in meters corresponding to the input geodetic locations. The array is the same size as the input
            array except the last dimension, which is still of size 3, is now reserved for X=0, Y=1, Z=2.

        '''

        if type(llh) is not np.ndarray:
            llh = np.array(llh, dtype=float)

        lat = np.radians(llh[..., 0:1])
        lon = np.radians(llh[..., 1:2])
        height_meters = llh[..., 2:3]

        cosphi = np.cos(lat)
        sinphi = np.sin(lat)
        coslam = np.cos(lon)
        sinlam = np.sin(lon)

        fm1 = (1.0 - self.F)
        fm12 = fm1 * fm1
        C = 1.0 / np.sqrt(cosphi * cosphi + fm12 * sinphi * sinphi)
        S = fm12 * C
        P = (self.A * C + height_meters) * cosphi

        result = np.empty_like(llh)
        result[..., 0:1] = P * coslam
        result[..., 1:2] = P * sinlam
        result[..., 2:3] = (self.A * S + height_meters) * sinphi
        return result

    # -----------------------------------------------------------------------------
    #               get_osculating_spheroid
    # -----------------------------------------------------------------------------
    def get_osculating_spheroid(self, geodeticlatitude: float, geodeticlongitude: float):
        '''
        Calculates the sphere that best fits the ellipsoid surface (height = 0)
        at the geodetic location. This algorithm best fits the spheroid in the
        latitudinal direction (ie North South). A better choice for the future may be to use a look vector. The
        osculating sphere is passed to radiative transfer models to "best" represent the true earth at the
        region of interest using a spherical system.

        Notes
        -----
        The radius of curvature is given by:-

        ..  math::

            R = \\frac{ \\left[1 + \\left(\\frac{dy}{dx}\\right)^2\\right]^{\\frac{3}{2}}}{\\frac{d^2y}{dx^2}}

        and for an ellipse:

        ..  math::

            R = \\frac{1}{ab} \\left[ \\frac{a^2}{b^2}y^2 + \\frac{b^2}{a^2}x^2\\right]^{\\frac{3}{2}}

        Parameters
        ----------
        geodeticlatitude: float
            The geodetic latitude in degrees of the location where we need to fit the osculating sphere

        geodeticlongitude: float
            The geodetic longitude in degrees of the location where we need to fit the osculating sphere.

        Returns
        -------
        Tuple[ radius:float, offset:Tuple[float, float]]:
            Returns a two element tuple containing the radius of the best fit osculating spheroid in the first element
            and the offset of the origin of the osculating sphere from the origin of the true geoid in the second element.
        '''

        xunit = np.array([math.cos(math.radians(geodeticlongitude)), math.sin(math.radians(geodeticlongitude)), 0.0])
        yunit = np.array([0.0, 0.0, 1.0])
        location = self.xyz_from_llh((geodeticlatitude, geodeticlongitude, 0.0))

        y0 = location[2]												            # Get the vertical component of the point on the surface of the geoid
        x0 = np.sqrt(np.square(location[0]) + np.square(location[1]))			    # Get the horizontal component of the point on the surface of the geoid
        a = self.A
        b = self.A * (1.0 - self.F)
        a2 = a * a
        b2 = b * b
        a2y0 = a2 * y0
        b2x0 = b2 * x0
        r = 1.0 / (a * b) * ((a2y0 * y0 / b2 + b2x0 * x0 / a2) ** 1.5)			    # Get the radius of curvature at this location
        theta = np.arctan2(a2y0, b2x0)										        # Get the angle of the gradient of the vertical at the surface of the geoid
        dx = r * np.cos(theta)												        # Get the X offset of the center of curvature from the surface point.
        dy = r * np.sin(theta)												        # Get the Y offset of the center of curvature from the surface point.
        offset = location - dy * yunit - dx * xunit
        radius = r
        return (radius, offset)

    # -----------------------------------------------------------------------------
    #               xyz_west_south_up
    # -----------------------------------------------------------------------------
    def xyz_west_south_up(self, xyz=None, llh=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Returns the three right-handed, orthogonal unit vectors, `west`, `south` and `up` at the provided locations. The unit vectors
        are always given in xyz geocentric coordinates. The locations at which the unit vectors are required can be specified
        in either geocentric xyz coordinates using the keyword ``xyz`` or in geodetic latitude, longitude, height using keyword
        ``llh``.

        The algorithm allows all the locations to be specified using multi-dimensional arrays with the last dimension set
        to a size of 3 ali_elements. This allows the code to execute all location calculations in one pass of the code.

        It is not recommended, but it is possible to specify locations using both keywords in which case the user must ensure they are identically
        sized arrays and they have the geocentric and geodetic coordinates pertaining to the same exact locations

        Parameters
        ----------
        xyz: array[..., 3], optional
            The geocentric, xyz, locations in meters where the three unit vectors will be calculated. It can be multidimensional
            with the last dimension equal to 3, [0=X, 1=Y, 2=Z]. It is usually left undefined if keyword ``llh`` is used.

        llh: array[..., 3], optional
            An array of look unit vectors given in geodetic lat, lon, height (derees, degrees, meters). It can be
            multi-dimensional but the last dimension must be equal to 3, [0=lat, 1=long, 2=height]. It is usually left
            undefined if keyword ``xyz`` is used.

        Returns
        -------
        Tuple[ west:array[..., 3], south:array[..., 3], up:array[..., 3]]
            A three element tuple storing the west, south and up unit vectors expressed as dimensional numbers in geocentric,
            xyz coordinates. The array sizes are the same as the arrays used to define the input locations.

        Notes
        -----
        The directions of west and south are undefined at either pole. It is technically possible that NaN values could be returned
        at the pole, which happens if the horizontal vector at the pole is calculated to be exactly 0.0. In practice,
        this very rarely happens even when you put locations at the pole as there is usually sufficient numerical round off
        that keeps the horizontal component very small but not zero. But it can happen and its not checked for.

        See Also
        --------
        :meth:`xyz_north_east_down`
        '''
        location: np.ndarray = xyz
        geo: np.ndarray = llh

        if (location is None):
            assert (geo is not None)
            location = self.xyz_from_llh(geo)

        if (geo is None):
            assert (location is not None)
            geo = self.llh_from_xyz(location)

        lat = np.radians(geo[..., 0:1])
        coslat = np.cos(lat)
        sinlat = np.sin(lat)

        horiz = np.empty_like(location)                                             # vector in equatorial plane and observers meridian
        horiz[..., 0:1] = location[..., 0:1]                                        # Get the X coordinate of the geocentric location
        horiz[..., 1:2] = location[..., 1:2]                                        # Get the Y coordinate of the geocentric location
        horiz[..., 2:3] = 0.0
        horiz = _unitvector(horiz)
        vertical = np.zeros_like(location)                                           # unit vector in equatorial plane and observers meridian
        vertical[..., 2:3] = 1.0

        up = horiz * coslat + vertical * sinlat	                                    # get vertical  unit vector at location
        south = horiz * sinlat - vertical * coslat	                                # get due South unit vector at observer's location
        west = np.cross(south, up, axisa=-1, axisb=-1)                              # get due West  unit vector at observer's location
        return (west, south, up)

    # -----------------------------------------------------------------------------
    #               xyz_north_east_down
    # -----------------------------------------------------------------------------
    def xyz_north_east_down(self, xyz=None, llh=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Returns the three right-handed, orthogonal unit vectors, `north`, `east` and `down` at the provided locations.
        It is a small wrapper around :meth:`xyz_west_south_up` and all the parameters are described in that function
        with the returning variables being changed from (``west``, ``south``,  ``up``) to (``north``, ``east``, ``down``)
        '''

        west, south, up = self.xyz_west_south_up(xyz=xyz, llh=llh)
        return (-south, -west, -up)

    # -----------------------------------------------------------------------------
    #               _check_tangent_pt_convergence_old
    # -----------------------------------------------------------------------------
    def _check_tangent_pt_convergence_old(self, nominaltanptxyz: np.ndarray, lookxyz: np.ndarray, startdistance: float) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Internal routine used to diagnose the tangent point location code. This algorithm assumes we have
        a decent approximation o fthe tangent point, say within a few meters. It then steps a little bit in front and behind
        and looks at the variation of altitude. The tangent point should occur at the minimum altitude.

        The code tries to be efficient and fits a parabola to height returned from three points at (+l, 0, -l) from the
        nominal tangent point. It can fit the quadratic coefficients and getthe location of the quadratic peak.

        Unfortunately the code seems too clever for its own good. It actually made some quite bad estimates, moving by 2 or 3 meters when
        it probably should have stayed within 0.2 meters. The code needs to be debugged and fixed or re-written so it just
        steps through, brute-force style, looking for the mnimum height.

        '''

        locxyz = np.copy(nominaltanptxyz)                           # make a copy of the tangent point
        llh = self.llh_from_xyz(locxyz)                             # get the lat, lng, height of the tangent point
        l = startdistance                                           # assume we are within 128 meters
        h0 = llh[..., 2:3]                                          # get the height of the current tangent point
        for n in range(8):                                          # then for 8 iterations, half the search space each time
            loc1 = locxyz + lookxyz * l                             # Step to a point distance "l" after tangent point
            loc2 = locxyz - lookxyz * l                             # Step to a point distance "l" before tangent point
            llh1 = self.llh_from_xyz(loc1)                          # get the lat/lon/height at the first
            llh2 = self.llh_from_xyz(loc2)                          # and second location
            h1 = llh1[..., 2:3] - h0                                # And calculate the change in height
            h2 = llh2[..., 2:3] - h0                                # at both locations
            #                a = (h1 + h2)/(2* l *l)                # fit a parabola (al^2 + bl + c = deltaH),  c = 0 as deltaH = 0, when L = 0
            #                b = (h1 - h2)/(2 * l)                  # nice quick and dirty fit for this special config
            #                pk = -b/(2.0 * a)                      # and standard formula for the peak
            pk = (h2 - h1) / (h1 + h2) * (0.5 * l)                  # alternatively , just go straight to the answer
            locxyz = locxyz + pk * lookxyz                          # update the tangent point location
            llh = self.llh_from_xyz(locxyz)                         # and its lat, lng, height
            h0 = llh[..., 2:3]                                      # get the height of the new tangent point
            l = l / 2.0                                             # and go to 1/2 step size
            if l < 0.003:                                           # anything less than 3 mm is just noise in the system
                break
        delta = np.linalg.norm(locxyz - nominaltanptxyz, axis=-1)
        return (locxyz, delta)

    # -----------------------------------------------------------------------------
    #               _check_tangent_pt_convergence
    # -----------------------------------------------------------------------------
    def _check_tangent_pt_convergence(self, nominaltanptxyz: np.ndarray, lookxyz: np.ndarray, checkconvergence) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Internal routine used to diagnose the tangent point location code. This algrotithm assumes we have
        a decent approximation o fthe tangent point, say within a few meters. It then steps a little bit in front and behind
        and looks at the variation of altitude. The tangent point should occur at the minimum altitude.

        The code tries to be efficient and fits a parabola to height returned from three points at (+l, 0, -l) from the
        nominal tangent point. It can fit the quadratic coefficients and getthe location of the quadratic peak.

        Unfortunately the code seems too clever for its own good. It actually made some quite bad estimates, moving by 2 or 3 meters when
        it probably should have stayed within 0.2 meters. The code needs to be debugged and fixed or re-written so it just
        steps through, brute-force style, looking for the mnimum height.

        '''

        llh = self.llh_from_xyz(nominaltanptxyz)                    # get the lat, lng, height of the tangent point
        h0 = llh[..., 2:3]
        l = np.linspace(-1.0, 1.0, 1001)
        loc = lookxyz * l[:, np.newaxis]                            # Step to a point distance "l" after tangent point
        loc2 = nominaltanptxyz + loc
        llh2 = self.llh_from_xyz(loc2)                              # and second location
        h2 = (llh2[..., 2:3] - h0).flatten()                       # at both locations
        minidx = np.argmin(h2)
        coeffs = np.polyfit(l, h2, 2)
        a = coeffs[0]
        b = coeffs[1]
        pk = -b / (2.0 * a)
        print('Tangent point convergence seems to be {} meters from the true tangent point'.format(pk))
        return pk

    # -----------------------------------------------------------------------------
    #               _mask_out_bad_vectors
    # -----------------------------------------------------------------------------
    def _mask_out_bad_vectors(self, badflags: np.ndarray, vectors: np.ndarray):
        '''
        Masks out any bad vectors by setting the vector to NaN. Assume we have a multi-dimensional array of vectors: eg.  vectors [N,M,K...,L,3]
        where the last dimension is the XYZ or Lat,Long,height. We also have an multi-dimensional array of bad vectors
        of exactly the same size except the last dimension is 1: badflags[N,M,K...,L,3]

        The code uses np.nonzero to find the indices of the bad vectors as a set of tuple indices for each dimension.
        We then manually toggle the index of the last dimension.
        '''

        if np.any(badflags):
            badidx = np.nonzero(badflags)                               # Get a list of indices f where the vectors are bad. This is the same index as used for the X component
            vectors[badidx] = np.nan                                    # set the X component, index 0 to Nan

            for i, j in enumerate(badidx[-1]):
                badidx[-1][i] = 1                                       # Set the bad indices for the Y component.
            vectors[badidx] = np.nan                                    # Set the Y component to NaN

            for i, j in enumerate(badidx[-1]):
                badidx[-1][i] = 2                                       # Set the bad indices for the Z component.
            vectors[badidx] = np.nan                                    # Set the Z component to NaN
        return vectors

    # -----------------------------------------------------------------------------
    #               xyz_tangent_point_location
    # -----------------------------------------------------------------------------
    def xyz_tangent_point_location(self, observerxyz: np.ndarray, lookxyz: np.ndarray, checkconvergence=None, check_tangent_behind=False) -> np.ndarray:
        '''
        Calculates the xyz geocentric coordinates of the tangent point given a look vector away from an observers location
        towards the Earth.  The earth is modelled as an oblate spheroid.  All calculations assume straight ray paths. This
        is the second technique we developed to find the tangent point and seems to be the adhoc default. In practice it is
        neither substantially quicker or more accurate than the original technique. It seems to work well at sub-meter accuracy
        (typically around 0.1 meter accuracy), but beyond that more careful handling of numeric roundoff may be required.

        The code is implemented to process multiple look vectors and observer positions in one pass of the code using
        multi-dimensional arrays. However both the observer and look arrays must be the same shape as the same
        shape flattening factor is iteratively applied to both fields;  it fails if they have different shapes. The
        legacy version of this function can handle  broadcasting, eg a single observer and multiple lines of sight.

        Algorithm Details
        -----------------
        This technique performs an affine transformation (the ratio of distances
        is preserved under the transformation) that maps the oblate spheroid into a true sphere.
        The tangent point, using a transformed observer and look vector, is quickly found for a sphere using simple geometry.
        The transformed tangent point can then be shifted back to the oblate spheroid system.

        The technique is elegant but does have an achilles heal that we ideally need to use the oblate spheroid that goes
        through the tangent point rather than the one that follows the Earth's surface. This results in an
        iterative approach where the required oblate spheroid shape is estimated at the end of each iteration and the
        fit repeated. It seems to work okay but there are currently no tests for convergence and there is no theoretical basis
        showing that the new choice of flattening factor between iterations is correct. For example, its not even clear to me that the
        oblate spheroid for 40 km is parallel to the oblate spheroid at the ground and that the vertical vector at the ground
        is perpendicular to the oblate spheroid at 40 km. All of which play into the definition of tangent point.

        Limitations
        -----------
        The lookxyz and observerxyz vector must not be parallel.  If they are, then there is a whole family
        of solutions (probably the great circle) around the oblate spheroid. This condition is not tested for
        and may produce very stange results.

        Parameters
        ----------
        observerxyz: array[...,3]
            An array of observer positions given in geocentric xyz coordinates (meters). It can be multidimensional with the last dimension equal to 3, [0=X, 1=Y, 2=Z].
            It should have the same shape as variable ``lookxyz``

        lookxyz: array[...,3]
            An array of look unit vectors given in geocentric xyz coordinates (dimensionless). The look vectors must be away from
            the observer towards the Earth. It can be multidimensional with the last dimension equal to 3, [0=X, 1=Y, 2=Z]. Its
            shape must be the same as ``observerxyz``.

        checkconvergence: float or None
            Diagnostic internal flag used to provide information about the true tangent point assuming its close to this solution.
            Its reserved for internal usage.

        Returns
        -------
        array[...,3]
            The location of the tangent point in geocentric, xyz coordinates (meters). The final array size is determined
            from the broadcasting of ``observerxyz`` with ``lookxyz``, but the last dimension will be 3. Locations that failed
            to find a tangent point are set to Nan. This typically occurs if the look vector is not towards the Earths limb.

        See Also
        --------
        :meth:`xyz_tangent_point_location_legacy`
        '''

        F = self.F
        obsTransf = np.array(observerxyz)
        ellTransf = np.array(lookxyz)

        for n in range(5):
            transformFactor = 1.0 / (1.0 - F)
            obsTransf[..., 2:3] = observerxyz[..., 2:3] * transformFactor       # apply the transform to the observer
            ellTransf[..., 2:3] = lookxyz[..., 2:3] * transformFactor           # and the look vector
            ellt = _unitvector(ellTransf)                                       # make the look vector into a unit vector
            s = -_dotproduct(ellt, obsTransf)                                   # l = Robs.cos(theta)
            locxyz = obsTransf + s * ellt                                       # get the tangent point in the transformed system
            locxyz[..., 2:3] /= transformFactor                                 # then map back to the regular system
            llh = self.llh_from_xyz(locxyz)
            F = self.F * (self.A / (self.A + llh[..., 2:3]))                    # We have to use a flattening ratio that corresponds to the oblate spheroid at the tangent altitude. I'm not convinced this is completely correct

        if check_tangent_behind:
            bad = s < 0.0
            locxyz = self._mask_out_bad_vectors(bad, locxyz)
        if (checkconvergence is not None):
            answer = self._check_tangent_pt_convergence(locxyz, lookxyz, checkconvergence)
        else:
            answer = locxyz
        return answer

    # -----------------------------------------------------------------------------
    #               xyz_tangent_point_location_legacy
    # -----------------------------------------------------------------------------
    def xyz_tangent_point_location_legacy(self, observerxyz: np.ndarray, lookxyz: np.ndarray, checkconvergence=None, do_high_precision: bool=True):
        '''
        Calculates the xyz geocentric coordinates of the tangent point given a look vector away from an observers location
        towards the Earth.  The earth is modelled as an oblate spheroid.  All calculations assume straight ray paths. This
        is the first technique we developed to find the tangent point. It seems to work well at sub-meter accuracy
        (typically around 0.1 meter accuracy), but beyond that more careful handling of numeric roundoff may be required.

        The code is implemented to process multiple look vectors and observer positions in one pass of the code using
        multi-dimensional arrays and broadcasting.

        Algorithm Details
        -----------------
        Basic idea is to transform from the standard (primed) geocentric coordinate system
        to another system such that the look vector defines the X direction, the
        spacecraft position is a second vector that defines the plane that includes
        the straight line ray and the center of the earth.  The new Z axis is
        perpendicular to the plane.	Convenient because Z is exactly zero in the new
        coordinate frame.  Also convenient is the fact that the straight line ray
        is given by the equation y = yr  where yr is the (constant) Y coordinate
        of the spacecraft in the new reference frame.

        First find the equation of the oblate spheroid projected onto this plane.
        Second define the tangent point as the place on the oblate spheroid that
        is parallel to the ray direction (which is aligned in the X direction)

        Transform the tangent point back to the normal coordinate system.

        Limitations
        -----------
        The lookxyz and observerxyz vector must not be parallel.  If they are, then there is a whole family
        of solutions (probably the great circle) around the oblate spheroid. This condition is not tested for
        and may produce very stange results.

        Parameters
        ----------
        observerxyz: array[...,3]
            An array of observer positions given in geocentric xyz coordinates (meters). It can be multidimensional with the last dimension equal to 3, [0=X, 1=Y, 2=Z]

        lookxyz: array[...,3]
            An array of look unit vectors given in geocentric xyz coordinates (dimensionless). The look vectors must be away from
            the observer towards the Earth. It can be multidimensional with the last dimension equal to 3, [0=X, 1=Y, 2=Z]. Its
            shape must be the same as ``observerxyz`` or at least be broadcast compatible.

        checkconvergence: float or None
            Diagnostic internal flag used to provide information about the true tangent point assuming its close to this solution.
            Its reserved for internal usage.

        do_high_precision: bool
            A flag that requests that the calculation is iterated a second time to give it high precision (sub-millimeter).
            By default this option is enabled (True). It can be disabled for low precision calculations that want to go a little faster
            and are happy with height precision around ~0.5 m.

        Returns
        -------
        array[...,3]
            The location of the tangent point in geocentric, xyz coordinates (meters). The final array size is determined
            from the broadcasting of ``observerxyz`` with ``lookxyz``, but the last dimension will be 3. Locations that failed
            to find a tangent point are set to Nan. This typically occurs if the look vector is not towards the Earths limb.

        See Also
        --------
        :meth:`xyz_tangent_point_location`

        '''

        a = self.semi_major_axis
        c = self.semi_minor_axis                                        # semi minor axis (for earth's Z axis)
        west, south, up = self.xyz_west_south_up(xyz=observerxyz)       # get the local "UP" vector at the observer

        xunit = _unitvector(lookxyz)					                # define X unit vector as parallel to the look vector
        zunit = _unitvector(np.cross(lookxyz, observerxyz))		        # define Z unit vector as the ray tracing plane as perpendicular to the plane defined by the ray and the position vector.
        yunit = _unitvector(np.cross(zunit, xunit))			            # define Y unit vector so it forms a right angled system.

        w11 = xunit[..., 0:1]						                    # get the transformation numbers from the unit vectors
        w12 = yunit[..., 0:1]						                    # we have to do the transpose to get the proper array (primes are Earth coords, unprimed are plane coords)
        w21 = xunit[..., 1:2]						                    # x' = w11.x + w12.y + w13.z
        w22 = yunit[..., 1:2]						                    # y' = w21.x + w22.y + w23.z
        w31 = xunit[..., 2:3]						                    # z' = w31.x + w32.y + w33.z
        w32 = yunit[..., 2:3]						                    # for our calculations z is exactly zero.

        a2 = a * a
        c2 = c * c

        A = (np.square(w11) + np.square(w21)) / a2 + np.square(w31) / c2		# Calculate coefficients of the equation of oblate spheroid
        B = (np.square(w12) + np.square(w22)) / a2 + np.square(w32) / c2        # projected onto the straight line ray
        C = 2.0 * ((w11 * w12 + w21 * w22) / a2 + (w31 * w32) / c2)				# plane

        factor = 4.0 * A * B - C * C
        Ty = _dotproduct(observerxyz, yunit)						            # y coordinate of observer
        Tx = -C / A * np.sqrt(A / factor)                                       # x coordinate
        tanpt_xyz = (Tx * xunit) + (Ty * yunit)		                            # Now get Geocentric location in original coordinate system

        if (do_high_precision):                                                  # Apply a second iteration
            llh = self.llh_from_xyz(tanpt_xyz)                                  # That calculates the height using the ground surface
            h = llh[..., 2:3]                                                   # extracts the height
            a2 = (a + h) * (a + h)                                              # now re-adjust the semi-major and semi-minor axes so the "ground" is now at the nominal height
            c2 = (c + h) * (c + h)                                              # and repeat the calculation. This should iterate very rapidly for all atmospheric altitudes.

            A = (np.square(w11) + np.square(w21)) / a2 + np.square(w31) / c2    # Calculate coefficients of the equation of oblate spheroid
            B = (np.square(w12) + np.square(w22)) / a2 + np.square(w32) / c2    # projected onto the straight line ray
            C = 2.0 * ((w11 * w12 + w21 * w22) / a2 + (w31 * w32) / c2)         # plane

            factor = 4.0 * A * B - C * C
            Tx = -C / A * np.sqrt(A / factor)                                   # x coordinate
            tanpt_xyz = (Tx * xunit) + (Ty * yunit)                             # Now get Geocentric location in original coordinate system

        bad = _dotproduct(lookxyz, observerxyz) > 0.0                           # Look for any point looking away from Earth
        tanpt_xyz = self._mask_out_bad_vectors(bad, tanpt_xyz)                  # and set them to NaN if they are.

        if (checkconvergence is not None):                                      # check convergence if requested. This is internal and not supported.
            answer = self._check_tangent_pt_convergence(tanpt_xyz, lookxyz, checkconvergence)
        else:
            answer = tanpt_xyz
        return answer

    # -----------------------------------------------------------------------------
    #               xyz_altitude_intercepts
    # -----------------------------------------------------------------------------
    def xyz_altitude_intercepts(self, Hnd: np.ndarray, observerxyznd: np.ndarray, lookxyznd) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Calculate the entrance and exit points of a look vector with a given shell at height H. The routine allows the
        user to pass in a multi-dimensional arrays of look vector, observer position and shell height. The calculation is
        only meaningful in the look vector is towards the Earth and the required shell height is at or above the
        tangent point, i.e. the shell height must lie between the tangent point and the observer.

        The entrance and exit points are approximately symmetric about the tangent point.  The entrance point is located
        towards the observer on the near-side of the tangent point. The exit point is on the far side, furthest from the
        tangent point.

        The algorithm is internally implemented to iterate over the arrays and process the intercepts one case at a time.
        It currently uses ``scipy.optimize.zbrentq`` to find the shell heights which is not conducive to array based operations.
        Consequently, the method, simply due to its array looping implementation, may be substantially slower than the other methods in this
        object. This should not be a significant issue for hundreds or thousands of rays but will be an issue for millions of rays.

        Parameters
        ----------
        Hnd: np.ndarray
            The height in meters of the target shell intercepts. This is typically given as an array with the same shape as ``observerxynd`` and ``lookxynd``
            except the last dimension is 1 rather than 3. It is also allowed to be a scalar or a simple sequence for
            situations where that is compatible with the other parameters

        observerxyznd: np.ndarray
            A multi-dimensional array that provides the geocentric, xyz location of the observer in meters. The array can be
            any shape but the last dimension must be 3. It is recommended that this array has the same shape as ``lookxyznd``.
            The observer is typically a spacecraft or high altitude platform.

        lookxyznd: np.ndarray
            A multi-dimensional array that provides the geocentric, xyz, look vector of each ray as a dimensionless
            unit vector. The rays are treated as being away from the observer, towards the Earth, which is opposite to the actual direction
            of light propagation. The array can be any shape but the last dimension must be 3. It is recommended that this
            array has the same shape as ``observerxyznd``.

        Returns
        -------
        Tuple[entrypoint: np.ndarray, exitpoint: np.ndarray, ok: np.ndarray]
            Three element tuple containing (i) the location of the entry point of the ray into the shell, (ii) the location of the
            exit point of the ray from the shell and (iii) a boolean status flag array. All variables are the same shape
            as the input arrays.
        '''

        if type(Hnd) is not np.ndarray:
            if isinstance(Hnd, numbers.Number):
                Hnd = np.array([Hnd])
            else:
                Hnd = np.array(Hnd)

        N = lookxyznd.size // 3
        ok = (lookxyznd.size == observerxyznd.size) and (N == Hnd.size)
        if not ok:
            raise ValueError('We can only support arrays where they all have the same number of observations [{}]'.format(N))

        entrypointanswers = np.zeros_like(lookxyznd)
        exitpointanswers = np.zeros_like(lookxyznd)
        dims = list(lookxyznd.shape)
        dims[-1] = 1
        status = np.ones(dims, dtype=bool)

        lookiter = np.nditer(lookxyznd)
        obsiter = np.nditer(observerxyznd)
        hiter = np.nditer(Hnd)
        entrypoint_iter = np.nditer(entrypointanswers, op_flags=['readwrite'])
        exitpoint_iter = np.nditer(exitpointanswers, op_flags=['readwrite'])
        status_iter = np.nditer(status, op_flags=['readwrite'])

        lookxyz = np.zeros([3])
        observerxyz = np.zeros([3])
        for iloop in range(N):                              # Loop over the number of lines of sight
            for i in range(3):
                lookxyz[i] = float(lookiter.value)
                observerxyz[i] = float(obsiter.value)
                lookiter.iternext()
                obsiter.iternext()
            H = float(hiter.value)
            hiter.iternext()

            # --- now process on observerxyz, lookxyz

            tanptxyz = self.xyz_tangent_point_location(observerxyz, lookxyz, check_tangent_behind=False)    # Get the tangent point in xyz
            tanptllh = self.llh_from_xyz(tanptxyz)                                                          # and in lat, lon, alt
            Ht = tanptllh[2]														                        # get the tangent height

            entrypoint = np.zeros([3])
            exitpoint = np.zeros([3])
            ok = (H > Ht)																                    # and make sure there is a valid solution
            if (ok):																		                # if there is then
                if (self.F == 0.0):															                # If we have a pure sphere then the solution
                    lmax = np.sqrt((2.0 * self.A + H + Ht) * (H - Ht))								        # is considerably simplified
                    d = lmax * lookxyz
                    entrypoint = tanptxyz - d
                    exitpoint = tanptxyz + d
                else:
                    evaluator = _HeightOffsetEvaluator(self, tanptxyz, lookxyz, H)		                    # Create the height offset evaluator object
                    Rmin = self.semi_minor_axis													            # Get the Semi Minor axis
                    Rmax = self.semi_major_axis													            # Get the Semi major axis
                    lmin = 0.9 * np.sqrt((H - Ht) * (H + Ht + 2 * Rmin))							        # Minimum distance from tangent point to shell height
                    lmax = 1.1 * np.sqrt((H - Ht) * (H + Ht + 2 * Rmax))							        # Maximum distance from tangent point to intersection

                    entrypoint, ok = evaluator.find_crossing_point(lmin, lmax)			                    # Find the shell from tangent point towards observer
                    if (ok):
                        exitpoint, ok = evaluator.find_crossing_point(-lmin, -lmax)		                    # Find the shell from tangent point away from observer

            for i in range(3):
                entrypoint_iter[0] = entrypoint[i]
                exitpoint_iter[0] = exitpoint[i]
                entrypoint_iter.iternext()
                exitpoint_iter.iternext()
            status_iter[0] = ok
            status_iter.iternext()

        return entrypointanswers, exitpointanswers, status

    # -----------------------------------------------------------------------------
    #               xyz_lookvector_from_tangent_altitude
    # -----------------------------------------------------------------------------
    def xyz_lookvector_from_tangent_altitude(self, h: float, observerxyz: np.ndarray, boresightazimuthxyz: Union[np.ndarray, float], convergence_meters=0.1):
        '''
        Finds the look vector towards a given tangent altitude at a given azimuth at the observer location.  The observer
        is typically a spacecraft position and should be (well) above the requested tangent point (e.g. more than 5 km).
        The algorithm makes a  first approximation using the osculating sphere which may fail if the requested tangent
        altitude is too close to the observers altitude.

        An iterative method is used to find the look vector required for the given tangent altitude. The algorithm uses
        oblate spheroid Earth and straight line geometry. Solution iterates until calculated tangent altitude is within
        ``convergence_meters``  meters of the requested height.

        Parameters
        ----------
        h: float
            The requested tangent height in meters above the surface of the Earth. This must be a scalar

        observerxyx: np.ndarray [3,]
            The geocentric, xyz, location of the observer in meters.

        boresightazimuthxyz: np.ndarray[3,] or float
            The horizontal unit vector at the observer's location that defines the azimuth of the required look vector.
            The parameter can also be specified as a floating point number that gives the compass bearing in degrees
            (0= North, 90= East etc.) of the desired unit vector.

        convergence_meters: float, optional
            Specifies the convergence. The returned solution will be within this height distance in meters from the target tangent point
            It may prove difficult to get consistency much better than 0.1 meters due to computational round off issues. default is 0.1 meters

        Returns
        -------
        np.ndarray[3,]
            The geocentric, xyz, dimensionless unit vector that looks from the observer toward the target tangent point at the
            given horizontal azimuth.
        '''

        obsllh = self.llh_from_xyz(observerxyz)                                     # Get the observers lat/lon, height

        if isinstance(boresightazimuthxyz, numbers.Number):
            north, east, down = self.xyz_north_east_down(llh=obsllh)                # Get the local geographic unit vectors
            azi = np.radians(boresightazimuthxyz)                                   # set the azimuth o fthe look vector
            boresightazimuthxyz = north * np.cos(azi) + east * np.sin(azi)          # and setup the horizontal look vector direction

        xunit = _unitvector(observerxyz)									        # Get unit vector to the observer
        yunit = boresightazimuthxyz - _dotproduct(boresightazimuthxyz, xunit) * xunit 	# Get unit vector perpendicular to observer and the observation plane.)

        re_geocentric, offset = self.get_osculating_spheroid(obsllh[0], obsllh[1])  # Get the osculating sphere suitable for this observer
        radius = np.linalg.norm(observerxyz - offset)                               # and get the radius of the observer, after adjusting for offset
        costheta = (re_geocentric + h) / radius									    # Get angle to tangent altitude on spherical Earth
        dh = np.nan
        ok = (costheta <= 1.0)
        if (not ok):
            raise RuntimeError('The target tangent altitude look vector must be below the observer')
        else:
            theta = np.arccos(costheta)
            sintheta = np.sin(theta)
            lookv = yunit * costheta - xunit * sintheta						            # and initialize the look vector
            ok = False
            numloops = 0
            while (not ok) and (numloops < 100):
                tanxyz = self.xyz_tangent_point_location(observerxyz, lookv)		    # calculate the tangent point
                tanllh = self.llh_from_xyz(tanxyz)                                      # and convert to lat,lon, height
                newh = tanllh[2]													    # get the geodetic height of this tangent point
                dh = (newh - h)													        # get the difference in height from guess to geodetic
                ok = (np.abs(dh) < convergence_meters)									# see if we have converged
                if not ok:																# we haven't
                    dtheta = 0.8 * dh / (radius * sintheta)							    # adjust the estimate of theta
                    theta += dtheta													    # get a new theta
                    costheta = np.cos(theta)
                    sintheta = np.sin(theta)
                    lookv = yunit * costheta - xunit * sintheta						    # and unpdate the unit vectors
                    numloops += 1
            if (not ok):
                raise RuntimeWarning('Stopped looking for the target tangent altitude after {} iterations. Convergence = {} meters'.format(numloops, dh))

        return lookv

    # -----------------------------------------------------------------------------
    #               xyz_observer_sun_look_for_given_sun_and_tangent_point
    # -----------------------------------------------------------------------------
    def xyz_observer_sun_look_for_given_sun_and_tangent_point(self, boresight_tanpt_llh: np.ndarray, omega_ocf: np.ndarray, lookazi_at_tanpt_deg: float=75.0, sza: float=None, ssa: float=None, saa: float=None, observer_height_m: float=500000.0):
        r"""
        Generates an observer who is looking along an instrument bore-sight at a given tangent point location. The Sun is chosen so
        it is at a given solar zenith angle and solar azimuth angle at that tangent point.  The observer is located at a given height.
        The location of the observer and sun are calculated that meet the specified requirements. A set of look vectors are generated
        using the direction angles :math:`\Omega` given in the optical control frame, centered around the boresight.

        Parameters
        ----------
        boresight_tanpt_llh: Tuple[lat, long , height]
            THe latitude, longitude and height of the tangent point that the boresight should look at
        omega_ocf: np.ndarray [N, 2]
            An array of N look directions specified using the :math:`\Omega_y` and :math:`\Omega_x` angles in the optical control
            frame. Briefly, :math:`\Omega_y` is upward elevation from the boresight and :math:`\Omega_x` is azimuth in an anti-clockwise
            direction from the bore sight, looking  away from the instrument towards port (left). All angles are given in degrees.
            The array is typically a 2-D array with the last dimension set to 2. :math:`\Omega_y` corresponds to the [...,  0] elements
            and :math:`\Omega_x` corresponds to thge [..., 1] elements. The variable can be a list or tuple, or anything else
            that can be easily coerced into an array. The variable can also be 1 dimensional (eg [N]) in which case it is coerced into
            a 2 Dimensional array with :math:`\Omega_x` equal to zero.
        lookazi_at_tanpt_deg : float
            The geographic azimuth at the tangent point of the look vector from the observer to the tangent point. It is specified in degrees 0=N, 90=E etc. Default is 75 degrees.
        sza : float
            Solar zenith angle in degrees at the desired tangent point. Default is None. And it must be set.
        saa : float
            Solar azimuth angle in degrees at the desired tangent point. Default is None. One of the keywords ``saa`` or ``ssa`` must be set to a valid number.
        ssa : float
            Solar scattering angle in degrees at the desired tangent point. Default is None. One of the keywords ``saa`` or ``ssa`` must be set to a valid number.
        observer_height : float
            The height of the observer above the surface in meters.  This must be higher than the requested tangent point. Default is 500,000.0

        Returns
        -------
        Dict[str, np.ndarray]
            A four element dictionary with fields ``observer``, ``sun``, ``los`` and ``icf_unitvectors``. The ``observer` element is the geocentric, xyz, location of the observer (in meters).
            The ``sun`` element is a geocentric, xyz, unit vector toward the sun. The ``los`` element is the array of look unit vectors away from the instrument
            corresponding to each of the look directions given in ``omega_ocf`. The look vectors are returned in the same shape as the
            look directions except the last dimension is now 3 instead of 2 and is the geocentric, xyz, unit vector of each look direction.
            The ``icf_unitvectors`` is a three element tuple containing three unit vectors, :math:`(\hat{x}_{icf}, \hat{y}_{icf}, \hat{z}_{icf}),
            associated with the :ref:`icf`. The three ICF unit vetors are return as geocentric, xyz, unit vectors.
        """

        if type(omega_ocf) is not np.ndarray:
            omega_ocf = np.array(omega_ocf)

        if omega_ocf.ndim > 1:
            omegay = np.radians(omega_ocf[..., 0, np.newaxis])
            omegax = np.radians(omega_ocf[..., 1, np.newaxis])
        else:
            omegay = np.radians(omega_ocf[:, np.newaxis])
            omegax = np.array(((0.0), (0.0), (0.0)))

        if (sza is None):
            raise ValueError('Keyword sza must be set to a valid value (0 to 180 degrees).')

        if (ssa is not None):
            if (saa is not None):
                raise ValueError('You cannot set both solar scattering angle (ssa) and solar azimuth angle (saa). Only choose one.')

            sza_r = np.radians(sza)
            ssa_r = np.radians(ssa)

            if np.cos(ssa_r) / np.sin(sza_r) > 1:
                raise ValueError('Unphysical geometry, the requested scattering angle is larger than the maximum value, {}'.format(np.degrees(np.arccos(np.sin(sza_r)))))
            else:
                saa = np.degrees(np.arccos(np.cos(ssa_r) / np.sin(sza_r))) + lookazi_at_tanpt_deg

        if (saa is None):
            raise ValueError('You must one keyword of either solar scattering angle (ssa) or solar azimuth angle (saa). Only choose one.')

        reftp = self.xyz_from_llh(boresight_tanpt_llh)                                                       # Get the location of the reference tangent point
        north, east, down = self.xyz_north_east_down(xyz=reftp)                                              # and get the geographic unit vectors at the tangent point.

        sunh = np.cos(np.deg2rad(saa)) * north + np.sin(np.deg2rad(saa)) * east                             # Find the sun unit vector
        sun = -np.cos(np.deg2rad(sza)) * down + np.sin(np.deg2rad(sza)) * sunh                              # at the tangent point location

        reflook = np.cos(np.deg2rad(lookazi_at_tanpt_deg)) * north + np.sin(np.deg2rad(lookazi_at_tanpt_deg)) * east        # look vector is tangential at tanegent point
        observer, exitpt, status = self.xyz_altitude_intercepts(observer_height_m, reftp, reflook)           # get the location of the observer
        west, south, up = self.xyz_west_south_up(xyz=observer)                                               # and get the geographic unit vectors at the tangent point.

        x_icf = reflook                                                                                     # get the x_icf coordinates of the bore-sight as the look vector to the desired tangent point
        y_icf = np.cross(up, x_icf)
        y_icf /= np.linalg.norm(y_icf, axis=-1)
        z_icf = np.cross(x_icf, y_icf)
        lookvectors = (np.cos(omegax) * x_icf + np.sin(omegax) * y_icf) * np.cos(omegay) + z_icf * np.sin(omegay)
        return {'observer': observer, 'sun': sun, 'los': lookvectors, 'icf_unitvectors': (x_icf, y_icf, z_icf)}
