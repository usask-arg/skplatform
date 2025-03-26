from typing import Tuple, Union
import numpy as np
import logging
# import sasktran as sk
import sktimeutils
from .rotationmatrix import RotationMatrix, UnitVectors
from ..geodesy import Geodetic


# -----------------------------------------------------------------------------
#           Class PlatformPointing
# -----------------------------------------------------------------------------
class PlatformPointing:

    """
    A class used to manage the pointing requirements of an **instrument** mounted on a **platform**. The primary objective
    is to

    * Mount an instrument on a platform, eg balloon, aircraft, satellite, ground.
    * Place and orient the platform at a location in space, on the ground or in the atmosphere.
    * Generate unit vectors at the front aperture of the instrument representing lines of sight and similar that can be
      used by the radiative transfer codes.

    The PlatformPointing class generates rotation matrices to transform between the different standard control frames as well as
    methods to mount the instrument on the platform, position the platform and rotate the platform.
    """
    # -----------------------------------------------------------------------------
    #           __init__
    # -----------------------------------------------------------------------------

    def __init__(self):
        self._LCF_Init_PCF: RotationMatrix = RotationMatrix(array=[[1, 0, 0], [0, -1, 0], [0, 0, -1]])                  # Initial matrix to get LCF properly aligned with PCF. LCFx is +PCFx, LCFy is -PCFy and LCFz is -PCFz
        self._PCF_Init_GCF: RotationMatrix = RotationMatrix(array=[[0, 1, 0], [1, 0, 0], [0, 0, -1]])                   # Initial matrix to get PCF properly aligned with GCF. PCF X is GCFy (due North), PCF y is GCFx (due East), PCF Z is -GCFz (upwards)
        self._ICF_to_LCF: RotationMatrix = RotationMatrix()                                                             # Instrument Rotation Matrix. Aligns the instrument control frame with the lab control frame,
        self._LCF_to_PCF: RotationMatrix = RotationMatrix()                                                             # Rotation Matrix from Instrument Control Frame to Platform Control Frame. Applies monutnig information
        self._PCF_to_GCF: RotationMatrix = RotationMatrix()                                                             # Rotation Matrix from Platform Control Frame to Geodetic Control Frame. E.G. Applies the Yaw, Pitch, Roll to orient the platform in local geodetic space.
        self._GCF_to_ECEF: RotationMatrix = RotationMatrix()                                                            # Rotation Matrix to convert local geodetic coordinates to geographic geocentric. ECEF
        self._utc: np.datetime64 = None                                                                                 #: Platform pointing also stores the UTC time as a convenience
        self._geo: Geodetic = Geodetic()                                                                                #: An instance of Geodetic object.
        self._geolocation: np.ndarray = np.zeros([3])
        self._local_west: np.ndarray = None
        self._local_south: np.ndarray = None
        self._local_up: np.ndarray = None
        self._nadirmode: bool = False

    # ------------------------------------------------------------------------------
    #           mount_instrument_on_platform
    # ------------------------------------------------------------------------------
    def mount_instrument_on_platform(self,
                                     azimuth_degrees: Union[np.ndarray, float],
                                     elevation_degrees: Union[np.ndarray, float],
                                     roll_degrees: Union[np.ndarray, float] = 0.0):
        """
        Mounts an instrument on a platform by setting up the azimuth, elevation and roll of the :ref:`icf`
        with respect to the :ref:`lcf`. Note that we try to use fairly standard definitions of azimuth, elevation
        and roll. The laboratory defines a local  X, Y, Z coordinate system. eg X is towards 0 degrees (or front) in a horizontal
        plane, Z is upwards and Y is in the horizontal plane and forms a right handed system towards the left/port.

        azimuth is left handed rotation around Z. It is clockwise rotation when viewed from above and follows a compass.
        elevation is left handed rotation around Y. It is positive upwards and is 0 in the horizontal plane
        roll is right handed rotation around X.

        Rotation are applied in the order Z, Y, X rotations in the order, azimuth first, elevation second and
        roll third,

        * azimuth is 0.0 when :math:`\\hat{x}_{ICF}` is pointing along the fuselage towards the nose,  i.e. it is parallel to :math:`\\hat{x}_{PCF}`. It is positive towards startboard (right hand).
        * elevation is 0.0 along the fuselage towards the nose when :math:`\\hat{x}_{ICF}` is parallel to :math:`\\hat{x}_{PCF}` and :math:`\\hat{z}_{ICF}` is pointing upwards. It is positive as :math:`\\hat{x}_{ICF}` tilts upwards.
        * roll is 0.0 when :math:`\\hat{y}_{ICF}` is parallel to the wings pointing to port-side (left) and is anti-parallel to :math:`\\hat{y}_{PCF}`. It is positive as :math:`\\hat{y}_{ICF}` tilts upwards and the aircraft leans to the right.

        The instrument control frame is initally configured so the instrument is mounted with azimuth, elevation and roll
        all equal to 0:  :math:`\\hat{x}_{ICF}` is parallel to :math:`\\hat{x}_{PCF}` and :math:`\\hat{z}_{ICF}` is
        anti-parallel to :math:`\\hat{z}_{PCF}`

        Parameters:
          azimuth_degrees : array[N] or float
             The azimuth of the instrument boresight in the platform control frame in degrees. Azimuth is zero along the aircraft
             fuselage towards the nose and increases towards the right/starboard
          elevation_degrees : array[N] or float
             The elevation of the instrument boresight in the platform control frame in degrees. Elevation is
             0.0 along the fuselage towards the nose and increases in the upward direction.
          roll_degrees : array[N] or float
             The right-handed roll in degrees of the instrument around the final resting point of the boresight, :math:`\\hat{x}_{ICF}`, in the platform control frame. The rotation
             is a right-handed roll around the final. Default is 0.0
        """

        self._ICF_to_LCF = RotationMatrix.from_yaw_pitch_roll(-np.radians(azimuth_degrees), -np.radians(elevation_degrees), np.radians(roll_degrees))

    # ------------------------------------------------------------------------------
    #               rotate_instrument_in_icf
    # ------------------------------------------------------------------------------
    def rotate_instrument_in_icf(self,
                                 azimuth_degrees: Union[np.ndarray, float],
                                 elevation_degrees: Union[np.ndarray, float],
                                 roll_degrees: Union[np.ndarray, float] = 0.0):
        """
        Sets the instrument rotation matrix so it rotates the instrument control frame, :math:`(\\hat{x}_{ICF}, \\hat{y}_{ICF}, \\hat{z}_{ICF})`,
        through the given azimuth, elevation and roll angles. This method is intended to simulate tilting mirrors and rotation stages
        attached to the instrument. The rotation is applied in the order, azimuth, elevation, roll.

        Parameters:
          azimuth_degrees : array[N] or float
             The azimuth of the instrument boresight in the platform control frame in degrees. Azimuth is left handed rotation
             around :math:`\\hat{z}_{ICF}`.
          elevation_degrees : array[N] or float
             The elevation of the instrument boresight in the platform control frame in degrees. Elevation is
             a left-handed rotation around :math:`\\hat{y}_{ICF}`.
          roll_degrees : array[N] or float
             The roll in degrees of the instrument. It is a right-handed rotation around :math:`\\hat{x}_{ICF}`. Default is 0.0
        """

        self._ICF_to_LCF = RotationMatrix.from_azimuth_elevation_roll(np.radians(azimuth_degrees), np.radians(elevation_degrees), np.radians(roll_degrees))

    # -----------------------------------------------------------------------------
    #           orient_platform_in_space
    # -----------------------------------------------------------------------------
    def orient_platform_in_space(self,
                                 yaw_degrees: Union[np.ndarray, float],
                                 pitch_degrees: Union[np.ndarray, float],
                                 roll_degrees: Union[np.ndarray, float]):
        """
        Sets the platforms yaw, pitch and roll to orient the platform in space at its location. The yaw, pitch
        and roll are applied in a geodetic system so angles are always about the local geodetic system regardless of
        the platform location and the values stay in effect until explicitly changed.

        This function is typically not required in simulators (although it can be set if you wish) as the platform pointing
        in these cases is overridden by the need to look at a "simulated" target location. The function is useful for
        attitude solutions from actual platforms where yaw, pitch and roll are a common method of specifying platform
        orientation

        This function will orient the entire platform in space. The platform control frame is initially set so the
        interpretation of yaw, pitch and roll are sensible, see :ref:`pcf`.

        Parameters:
          yaw_degrees : array[N] or float
             The yaw or geographic bearing of the platform control frame in degrees. This is a geographic bearing, N= 0, E= 90, S= 180, W = 270.
          pitch_degrees : array[N] or float
             The pitch of the platform control frame in degrees. Pitch is right-handed rotation around the starboard pointing :math:`\\hat{y}_{PCF}`.
             A positive pitch will tip the :math:`\\hat{x}_{PCF}` unit vector upwards. Negative values will tip it downwards.
          roll_degrees : array[N] or float
             The roll of the platform control frame in degrees. Roll is the right-handed rotation around :math:`\\hat{x}_{PCF}`
             after yaw and pitch have been applied. A pilot in a plane looking out of the nose along :math:`\\hat{x}_{PCF}`
             will roll to his right side, (right shoulder falls, left shoulder rises) assuming he is not flying upside down.

        """

        self._PCF_to_GCF = RotationMatrix.from_yaw_pitch_roll(np.radians(yaw_degrees), np.radians(pitch_degrees), np.radians(roll_degrees))

    # ------------------------------------------------------------------------------
    #           update_location
    # ------------------------------------------------------------------------------
    def set_platform_location(self,
                              xyzt: Tuple[Union[np.ndarray, float], Union[np.ndarray, float], Union[np.ndarray, float], Union[np.ndarray, np.datetime64]] = None,
                              latlonheightandt: Tuple[Union[np.ndarray, float], Union[np.ndarray, float], Union[np.ndarray, float], Union[np.ndarray, np.datetime64]] = None):
        """
        Sets the location and time of the platform to the given location and time. The methods supports setting just one location
        or an array of locations. The array option allows it to be used reasonably efficiently for instruments or simulations that need
        to process many (tens or hundreds of thousands of) locations.

        The location is used to configure the internal ``GCF_to_ECEF`` rotation matrix that transforms local geodetic
        look vectors into global ECEF X,Y,Z look vectors.

        The caller can choose to specify the location as either global geocentric ECEF coordinates using keyword option ``xyzt``
        or as geodetic latitude, longitude and altitude using keyword option ``latlonheight``. The caller must exlicitly set one
        of these options

        Parameters:
          xyzt: Tuple[float, float, float, np.datetime64] **or** Tuple[ array[N], array[N], array[N], array[N]]
              A 4 element tuple that specifies the platform location and time with  X, Y, Z  ECEF coordinates in meters and UTC

                 #. ECEF X geocentric coordinate in meters, either scalar or array[N]
                 #. ECEF Y geocentric corodinate in meters, either scalar or array[N]
                 #. ECEF Z geocentric coordinate in meters, either scalar or array[N]
                 #. UTC in any format that ``sktimeutils.ut_to_datetime64`` supports, either scalar or array[N]

          latlonheight: Tuple[float, float, float, np.datetime64] **or** Tuple[ array[N], array[N], array[N], array[N]]
              A 4 element tuple that specifies the platform location with  Geodetic coordinates X, Y, Z  ECEF geocentric location of the platform in meters.

                #. Geodetic latitude in degrees, either scalar or array[N]
                #. Geodetic longitude in degrees, either scalar or array[N]
                #. Geodetic altitude in meters above sea level, either scalar or array[N]
                #. UTC in any format that ``sktimeutils.ut_to_datetime64`` supports, either scalar or array[N]

        Returns:
           True if successful otherwise False.
         """

        ok = True

        if (xyzt is not None):
            t = xyzt[3]
            location = np.array([xyzt[0], xyzt[1], xyzt[2]])
        elif (latlonheightandt is not None):
            llh = np.array([latlonheightandt[0], latlonheightandt[1], latlonheightandt[2]])
            location = self._geo.xyz_from_llh(llh)
            t = latlonheightandt[3]
        else:
            raise ValueError("You must set either xyz or latlonheight")

        utc = sktimeutils.ut_to_datetime64(t)
        localwest, localsouth, localup = self._geo.xyz_west_south_up(xyz=location)

        isarraybased = (location.size // 3) > 1
        if isarraybased:
            self._GCF_to_ECEF = RotationMatrix.from_transform_to_destination_coordinates(-localwest, -localsouth, localup)
        else:
            self._GCF_to_ECEF = RotationMatrix.from_transform_to_destination_coordinates(-localwest, -localsouth, localup)

        # Keep these defined for the moment for backward compatibility. They only make good sense wrt to scalar calls
        self._geolocation = location[:]
        self._utc = utc
        self._local_west = localwest[:]
        self._local_south = localsouth[:]
        self._local_up = localup[:]
        # **** TODO **** Check this works with single points and multiple points
        return ok

    # -----------------------------------------------------------------------------
    #           location
    # -----------------------------------------------------------------------------
    def location(self):
        "returns the geocentric location of the platform. Only valid after a successfull call to set_platform_location"
        return self._geolocation

    # -----------------------------------------------------------------------------
    #           platform_utc
    # -----------------------------------------------------------------------------
    def utc(self):
        return self._utc

    # -----------------------------------------------------------------------------
    #           local_west
    # -----------------------------------------------------------------------------

    def local_west(self):
        """
        Returns the geocentric unit vector of west at the current location of the platform. Only valid after a successful call to set_platform_location
        """
        return self._local_west

    # -----------------------------------------------------------------------------
    #           local_south
    # -----------------------------------------------------------------------------

    def local_south(self):
        """
        Returns the geocentric unit vector of south at the current location of the platform. Only valid after a successful call to set_platform_location
        """
        return self._local_south

    # -----------------------------------------------------------------------------
    #           local_up
    # -----------------------------------------------------------------------------

    def local_up(self):
        """
        Returns the geocentric unit vector of up at the current location of the platform. Only valid after a successful call to set_platform_location
        """
        return self._local_up

    # -----------------------------------------------------------------------------
    #               reset_icf_rotation_matrices
    # -----------------------------------------------------------------------------
    def reset_icf_rotation_matrices(self):
        """
        Resets the IRM and ICF_to_PCF matrices to unity.
        """
        self._ICF_to_LCF = RotationMatrix(RotationMatrix.IUnit())
        self._LCF_to_PCF = RotationMatrix(RotationMatrix.IUnit())

    # -----------------------------------------------------------------------------
    #           force_pcf_rotation_matrix
    # -----------------------------------------------------------------------------
    def force_pcf_rotation_matrix(self, GEO: UnitVectors):
        """
        Defines the rotation matrix applied to the platform control frame so the :ref:`icf` unit vectors axes are
        aligned to the ECEF unit-vectors passed in after the full stack of rotation matrices are applied. Note this
        rotation is not applied to the Instrument Rotation Matrix (_ICF_to_LCF) as it is intended to force the primary
        boresight of the instrument control frame point towards a specific location.

        The problem is a linear algebra problem where we have G=RX, where G is the 3x3 Unitvector array passed in,
        X is the initial 3x3 unit vectors in the instrument control frame and R is the stack of rotation matrices. The
        matrix R can be expanded and we get a matrix equation of the form (A )(PCF_to_GCF) ( B )= G where we know
        everything except PCF_to_GCF. Thus we get  PCF_to_GCF = (A-1) G (B-1)

        Parameters:
          GEO : UnitVectors
             A 3x3 array of column unit vectors. These unit vectors specify the desired orientation of the instrument control frame
             vectors after rotation to the :ref:`ecef`.
        """
        B = RotationMatrix(self._PCF_Init_GCF.R  # B is the rotation matrics from instrument to Start of geodetic control frome
                           @ self._LCF_to_PCF.R
                           @ self._LCF_Init_PCF.R)

        A = self._GCF_to_ECEF                                                 # A is the rotation matrix from geodetic control frome to GEO

        self._PCF_to_GCF.from_rotation_matrix(A.RInv @ GEO.R @ B.RInv)       # get the platform "yaw,pitch roll " rotation matrix

    # -----------------------------------------------------------------------------
    #               convert_icf_to_ecef
    # -----------------------------------------------------------------------------
    def convert_icf_to_ecef(self, v: np.ndarray) -> np.ndarray:
        """
        Converts vectors expressed in the instrument coordinate frame to geographic geocentric ECEF vectors using the
        current rotation matrices.

        Parameters:
          v : np.ndarray
             An array (3xN) of N vectors expressed in the instrument control frame.

             * v[0,:] is the :math:`\\hat{x}_{ICF}` component of each vector,
             * v[1,:] is the :math:`\\hat{y}_{ICF}` component,
             * v[2,:] is the :math:`\\hat{z}_{ICF}` component.

        Returns:
          np.ndarray
            An array (3xN) of N vectors expressed in the geographic geocentric control frame,

            * v[0,:] is the :math:`\\hat{x}_{ECEF}` component,
            * v[1,:] is the :math:`\\hat{y}_{ECEF}` component,
            * v[2,:] is the :math:`\\hat{z}_{ECEF}` component.
        """
        R = self.get_icf_to_ecef_matrix()
        vnew = R.R @ v
        return vnew

    # ------------------------------------------------------------------------------
    #               convert_icf_to_ecef(self, v: np.ndarray)->np.ndarray:
    # ------------------------------------------------------------------------------

    def convert_icf_to_gcf(self, v: np.ndarray) -> np.ndarray:
        """
        Converts vectors expressed in the instrument coordinate frame to geodetic control frame vectors using the
        current rotation matrices.

        Parameters:
          v : np.ndarray
             An array (3xN) of N vectors expressed in the instrument control frame.

             * v[0,:] is the :math:`\\hat{x}_{ICF}` component of each vector,
             * v[1,:] is the :math:`\\hat{y}_{ICF}` component,
             * v[2,:] is the :math:`\\hat{z}_{ICF}` component.

        Returns:
          np.ndarray
            An array (3xN) of N vectors expressed in the geodetic control frame,

            * v[0,:] is the :math:`\\hat{x}_{GCF}` (West) component,
            * v[1,:] is the :math:`\\hat{y}_{GCF}` (South) component ,
            * v[2,:] is the :math:`\\hat{z}_{GCF}` (Up) component.
        """

        R = (self._PCF_to_GCF.R
             @ self._PCF_Init_GCF.R
             @ self._LCF_to_PCF.R
             @ self._LCF_Init_PCF.R
             @ self._ICF_to_LCF.R)
        vnew = R @ v
        return vnew

    # -----------------------------------------------------------------------------
    #           convert_gcf_to_ecef
    # -----------------------------------------------------------------------------

    def convert_gcf_to_ecef(self, v: np.ndarray) -> np.ndarray:
        """
        Converts vectors expressed in the geodetic control frame to geographic geocentric vectors using the
        current rotation matrices.

        Parameters:
          v : np.ndarray
             An array (3xN) of N vectors expressed in the geodetic control frame.

             * v[0,:] is the :math:`\\hat{x}_{GCF}` component of each vector,
             * v[1,:] is the :math:`\\hat{y}_{GCF}` component,
             * v[2,:] is the :math:`\\hat{z}_{GCF}` component.

        Returns:
          np.ndarray
             An array (3xN) of N vectors expressed in the geodetic control frame,

             * v[0,:] is the :math:`\\hat{x}_{ECEF}` component,
             * v[1,:] is the :math:`\\hat{y}_{ECEF}` component,
             * v[2,:] is the :math:`\\hat{z}_{ECEF}` component.
        """

        vnew = self._GCF_to_ECEF.R @ v
        return vnew

    # -----------------------------------------------------------------------------
    #           get_icf_to_ecef_matrix
    # -----------------------------------------------------------------------------
    def get_icf_to_ecef_matrix(self) -> RotationMatrix:
        """
        Returns the rotation matrix that converts vectors expressed in the instrument control frame to vectors in the geographic geocentric control frame
        If you are inspecting the source code for this method note that the
        rotation matrices are applied in reverse order, ie. right-most or last array in the list is the first rotation applied and left-most or first in the list is the last
        rotation/operation.

        Returns:
          np.ndarray
            An (3x3) rotation matrix. The array is used to transform vectors from the ICF instrument control frame to
            the ECEF geocentric control frame. The matrix should be applied as :math:`\\boldsymbol{V_{GEO}} = \\boldsymbol{R} @ \\boldsymbol{V_{ICF}}`
        """

        full_rotation = (self._GCF_to_ECEF.R                                                                            # 6) Convert the geodetic vectors to geocentric vectors.
                         @ self._PCF_to_GCF.R                                                                           # 5) Apply the spatial orientation of the platform in geodetic coordinates
                         @ self._PCF_Init_GCF.R                                                                         # 4) Rotate the platform control frame so its properly aligned with the geodetic corodinates
                         @ self._LCF_to_PCF.R                                                                           # 3) Rotate the instrument boresight so its properly positioned on the platform
                         @ self._LCF_Init_PCF.R                                                                         # 2) Flip lab control frame so it is in its initial orientation in the platform control frame
                         @ self._ICF_to_LCF.R)                                                                          # 1) Apply Instrument Control Frame to Lab Control Frame

        return RotationMatrix(array=full_rotation)
