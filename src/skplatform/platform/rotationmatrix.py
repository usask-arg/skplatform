from typing import Union, Tuple, List
import numbers
import numpy as np


# -----------------------------------------------------------------------------
#           class RotationMatrix:
# -----------------------------------------------------------------------------
class RotationMatrix:
    """
    A class that provides methods to calculate the rotation and transformation matrices used elsewhere in the package.
    These matrices are used to transform lines of sight and other vectors from one control frame to another.

        -  :py:meth:`from_transform_to_destination_coordinates`
        -  :py:meth:`from_rotate_one_setofaxes_to_another`
        -  :py:meth:`from_yaw_pitch_roll`
        -  :py:meth:`from_ypr_in_first_frame_to_ypr_in_second_frame`
        -  :py:meth:`transform_vectors`

    """

    # ------------------------------------------------------------------------------
    #           IUnit
    # ------------------------------------------------------------------------------
    @staticmethod
    def IUnit():
        """
        Returns the identity matrix, :math:`\\boldsymbol{IUnit}` as a 3x3 matrix
        """
        return np.array([[1.0, 0.0, 0.0],  #: The rotation matiopn. By default it is the Unit matrix *"i.e. no effect"
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]], dtype='float64')

    # -----------------------------------------------------------------------------
    #           __init__
    # -----------------------------------------------------------------------------

    def __init__(self, array: Union[List[List[float]], np.ndarray] = None):

        self._R: np.ndarray = None
        self._R: np.ndarray = RotationMatrix.IUnit() if array is None else np.array(array, dtype='float64')            #: The rotation matiopn. By default it is the Unit Identity matrix *"i.e. no effect"
        self._arraybased = False
        if (array is None):
            self._R = RotationMatrix.IUnit()
        else:
            self._R = array

    # -----------------------------------------------------------------------------
    #               __matmul__
    # -----------------------------------------------------------------------------
    def __matmul__(self, other):

        answer = RotationMatrix()
        if type(other) is RotationMatrix:
            answer._R = self._R @ other._R
        else:
            answer._R = self._R @ other
        return answer

    # -----------------------------------------------------------------------------
    #           property R
    # -----------------------------------------------------------------------------
    @property
    def R(self):
        return self._R

    @property
    def RInv(self):
        n = self._R.ndim
        if (n == 2):
            rinv = np.transpose(self._R, axes=(n - 1, n - 2))
        else:
            rinv = np.transpose(self._R, axes=(0, n - 1, n - 2))
        return rinv

    # -----------------------------------------------------------------------------
    #               _numpoints_in_mustbearrayargs
    # -----------------------------------------------------------------------------
    @staticmethod
    def _numpoints_in_mustbearrayargs(arguments: Tuple[Union[np.ndarray, float], ...]):
        '''
        Given a list of 1-D arrays and/or scalar numbers find the size of the 1-D arrays. All the arrays in the list
        should have the same size and must not be 0.
        '''

        N = 1
        arraybased = False
        for arg in arguments:
            if arg.ndim > 1:
                arraybased = True
                n = arg.shape[0]
                if (n != 1):
                    if (N == 1):
                        N = n
                    assert (n == N)

        return N, arraybased

    # -----------------------------------------------------------------------------
    #               _numpoints_in_args
    # -----------------------------------------------------------------------------
    @staticmethod
    def _numpoints_in_args(arguments: Tuple[Union[np.ndarray, float], ...]):
        '''
        Given a list of 1-D arrays and/or scalar numbers find the size of the 1-D arrays. All the arrays in the list
        should have the same size and must not be 0.
        '''

        N = 1
        arraybased = False
        for arg in arguments:
            if type(arg) is np.ndarray:
                assert (arg.ndim == 1)
                arraybased = True
                n = arg.size
                if (n != 1):
                    assert (n > 0)
                    if (N == 1):
                        N = n

                    assert (n == N)
            else:
                assert (isinstance(arg, numbers.Number))

        return N, arraybased

    # -----------------------------------------------------------------------------
    #           from_rotation_matrix
    # -----------------------------------------------------------------------------
    def from_rotation_matrix(self, other: np.ndarray):
        """
        Sets the rotation from another 3x3 matrix
        """
        self._R = other

    # -----------------------------------------------------------------------------
    #           from_transform_to_destination
    # -----------------------------------------------------------------------------
    @staticmethod
    def from_transform_to_destination_coordinates(ax: np.ndarray, ay: np.ndarray, az: np.ndarray):
        """
        Sets the rotation/transformation vector :math:`\\boldsymbol{R}` so it converts a vector, :math:`\\vec{v}`, specified in the system, :math:`\\boldsymbol{A}` (
        :math:`\\hat{x}_a` , :math:`\\hat{y}_a` , :math:`\\hat{z}_a` ),  to the same vector but specified in the destination system, :math:`\\boldsymbol{B}`
        (:math:`\\hat{x}_b` , :math:`\\hat{y}_b` , :math:`\\hat{z}_b`). This is implemented as a dot product.

        .. math::
               \\begin{split}
               \\vec{V}_A &= a_x\\hat{x}_a + a_y\\hat{y}_a + a_z\\hat{z}_a  \\\\
               \\vec{V}_B &= b_x\\hat{x}_b + b_y\\hat{y}_b + b_z\\hat{z}_b  \\\\
               \\end{split}

        where

        .. math::
               \\begin{split}
               \\hat{x}_a &=& a_{11}\\hat{x}_b +& a_{12}\\hat{y}_b +& a_{13}\\hat{z}_{b} \\\\
               \\hat{y}_a &=& a_{21}\\hat{x}_b +& a_{22}\\hat{y}_b +& a_{23}\\hat{z}_{b} \\\\
               \\hat{z}_a &=& a_{31}\\hat{x}_b +& a_{32}\\hat{y}_b +& a_{33}\\hat{z}_{b} \\\\
               \\end{split}

        then

        .. math::
                \\begin{pmatrix}
                Bx \\\\
                By \\\\
                Bz \\\\
                \\end{pmatrix}
                =
                \\begin{pmatrix}
                a_{11} & a_{21} & a_{31} \\\\
                a_{12} & a_{22} & a_{32} \\\\
                a_{13} & a_{23} & a_{33} \\\\
                \\end{pmatrix}
                \\begin{pmatrix}
                Ax  \\\\
                Ay  \\\\
                Az  \\\\
                \\end{pmatrix}

        Note that the 3x3 matrix is the *transpose* of the unit vector equations

        Parameters:
          ax
            A 3 element array expressing the :math:`\\hat{ax}` unit vector of frame :math:`\\boldsymbol{A}` in terms of the
            unit vectors in :math:`\\boldsymbol{B}`
          ay
            A 3 element array expressing the :math:`\\hat{ay}` unit vector of frame :math:`\\boldsymbol{A}` in terms of the
            unit vectors in :math:`\\boldsymbol{B}`
          az
            A 3 element array expressing the :math:`\\hat{az}` unit vector of frame :math:`\\boldsymbol{A}` in terms of the
            unit vectors in :math:`\\boldsymbol{B}`.

        """
        N, isarraybased = RotationMatrix._numpoints_in_mustbearrayargs((ax, ay, ax))
        if (isarraybased):
            R = np.empty([N, 3, 3], dtype='float64')
            R[:, :, 0] = ax
            R[:, :, 1] = ay
            R[:, :, 2] = az
        else:
            R = np.empty([3, 3], dtype='float64')
            R[:, 0] = ax
            R[:, 1] = ay
            R[:, 2] = az
        return RotationMatrix(R)

    # -----------------------------------------------------------------------------
    #           from_rotate_one_setofaxes_to_another
    # -----------------------------------------------------------------------------
    @staticmethod
    def from_rotate_one_setofaxes_to_another(A: np.ndarray, B: np.ndarray) -> 'RotationMatrix':
        """
        Generates the rotation matrix that rotates vectors specified in reference frame :math:`\\boldsymbol{A}` to vectors
        specified in reference frame :math:`\\boldsymbol{B}`. IUnit.E. the 3 unit vectors of frame :math:`\\boldsymbol{A}`  will be
        rotated on to the 3 unit vectors of frame  :math:`\\boldsymbol{B}`. The unit vectors for system :math:`\\boldsymbol{A}`, :math:`\\begin{bmatrix} \\hat{x}_A & \\hat{y}_A & \\hat{z}_A \\end{bmatrix}`
        and system :math:`\\boldsymbol{B}`, :math:`\\begin{bmatrix} \\hat{x}_B & \\hat{y}_B & \\hat{z}_B \\end{bmatrix}` are specified in a common
        system :math:`\\begin{bmatrix} \\hat{x'} & \\hat{y'} & \\hat{z'} \\end{bmatrix}` as 3x3 matrices.

        .. math::

                A= \\begin{pmatrix}
                   \\hat{x}_{A}.\\hat{x}' & \\hat{y}_{A}.\\hat{x}' & \\hat{z}_{A}.\\hat{x}' \\\\
                   \\hat{x}_{A}.\\hat{y}' & \\hat{y}_{A}.\\hat{y}' & \\hat{z}_{A}.\\hat{y}' \\\\
                   \\hat{x}_{A}.\\hat{z}' & \\hat{y}_{A}.\\hat{z}' & \\hat{z}_{A}.\\hat{z}' \\\\
                   \\end{pmatrix}

        .. math::

                B = \\begin{pmatrix}
                \\hat{x}_{B}.\\hat{x}' & \\hat{y}_{B}.\\hat{x}' & \\hat{z}_{B}.\\hat{x}' \\\\
                \\hat{x}_{B}.\\hat{y}' & \\hat{y}_{B}.\\hat{y}' & \\hat{z}_{B}.\\hat{y}' \\\\
                \\hat{x}_{B}.\\hat{z}' & \\hat{y}_{B}.\\hat{z}' & \\hat{z}_{B}.\\hat{z}' \\\\
                \\end{pmatrix}

        The rotation matrix :math:`\\boldsymbol{R}` is the operation that transforms the three unit vectors of :math:`\\boldsymbol{A}` in to
        the three unit vectors of :math:`\\boldsymbol{B}`

        .. math::

                \\boldsymbol{RA} &=  \\boldsymbol{B}\\\\
                \\boldsymbol{R}  &=  \\boldsymbol{BA^{-1}}\\\\

        Parameters:
          A
            A 3x3 matrix with the 3 unit vectors of control frame :math:`\\boldsymbol{A}` :math:`\\begin{bmatrix} \\hat{x}_A & \\hat{y}_A & \\hat{z}_A \\end{bmatrix}`
            expressed in a common coordinate sytem, :math:`\\hat{x}'`, :math:`\\hat{y}'`, :math:`\\hat{z}'`.

          B
            A 3x3 matrix with the 3 unit vectors of control frame :math:`\\boldsymbol{B}` :math:`\\begin{bmatrix} \\hat{x}_B & \\hat{y}_B & \\hat{z}_B \\end{bmatrix}`
            expressed in a common coordinate sytem, :math:`\\hat{x}'`, :math:`\\hat{y}'`, :math:`\\hat{z}'`.

        **Examples**:

        Rotate Platform and Instrument to look at point on Earth
            A typical usage example is to calculate the rotation matrix for a platform so an instrument on that platform
            will look at a point on the Earth. The solution is to align the platform control frame with the local geodetic
            system (west, south, up) and then express the instrument control frame, :math:`\\boldsymbol{A}` , as geographic
            geocentric unit vectors via the west,south, up unit vectors.  The :math:`\\boldsymbol{B}` unit vectors are taken
            as the look unit vector from the platform to the desired point on Earth, :math:`\\hat{x}_B`, while
            :math:`\\hat{y}_B` and :math:`\\hat{z}_B` are chosen to meet observation requirements, eg horizontal slit etc.
            In this example both :math:`\\boldsymbol{A}` and :math:`\\boldsymbol{B}` can be readily expressed in geographic
            geocentric coordinates using other parts of the software package.

        Rotate Instrument Control Frame to Platform Control Frame
            For this example we want to rotate vectors specified in the instrument control frame so they are specified
            in the platform control frame. For example the CATS instrument was located on the backside of the gondola
            during the 2018 campaign and required the line of sight to be rotated so it looked backwards. The solution
            is specify the instrument control frame unit vectors in platform control frame coordinates and recognise that the
            common coordinates are the platform control frame unit vectors. In this case the instrument control frame is system
            :math:`\\boldsymbol{A}` and the platform control frame is system :math:`\\boldsymbol{B}`. Because system :math:`\\boldsymbol{B}`
            is the common system then matrix :math:`\\boldsymbol{B }` is simply the unit vector :math:`\\boldsymbol{IUnit}`.
            Matrix :math:`\\boldsymbol{A}` must be calculated from instrument and platform configuration details but
            then :math:`\\boldsymbol{R}` is readily derived.


        """

        R = B @ np.linalg.inv(A)       # Note that you can also use the transpose of A, as it should be the inverse for proper rotations
        return RotationMatrix(R)

    # -----------------------------------------------------------------------------
    #           from_yaw_pitch_roll
    # -----------------------------------------------------------------------------
    @staticmethod
    def from_yaw_pitch_roll(yaw: Union[np.ndarray, float],
                            pitch: Union[np.ndarray, float],
                            roll: Union[np.ndarray, float]) -> 'RotationMatrix':
        """
        Creates a rotation matrix that applies yaw, pitch and roll in that order. This method is intended for aircraft
        system where :math:`\\hat{x}` is looking forward looking, :math:`\\hat{z}` is downwards and :math:`{y}` is towards starboard.
        All angles are made in a right-handed sense around the axis of rotation. Typically, positive pitch is turning "*upwards*",
        positive yaw is "*turning towards the right*" and positive roll is "*leaning to the right*".

        The matrix is created in the rotation order Yaw, Pitch, Roll (`rzyx`). The rotations are intrinsic, meaning the axes of
        rotations are also rotated. Intrinsic rotations require the rotations to be applied in reverse order (ie roll, pitch, yaw)
        but then the matrix operations are applied from right to left so we get the order yaw, @ pitch @ roll in the actual code to
        represent (yaw, pitch, roll).

        Parameters:
          yaw
            The yaw angle in radians. This rotation is applied first to the z axis.
          pitch
            The pitch angle in radians. This rotation is applied second to the rotated y axis.
          roll
            The roll angle in radians. This rotation is applied last to the rotated x axis.
        """

        N, isarraybased = RotationMatrix._numpoints_in_args((yaw, pitch, roll))
        yawMatrix = np.empty([N, 3, 3], dtype='float64')
        pitchMatrix = np.empty([N, 3, 3], dtype='float64')
        rollMatrix = np.empty([N, 3, 3], dtype='float64')

        yawMatrix[:, 0, 0] = np.cos(yaw)                                                                                # Yaw is rotation around the Z axis
        yawMatrix[:, 0, 1] = -np.sin(yaw)
        yawMatrix[:, 0, 2] = 0
        yawMatrix[:, 1, 0] = np.sin(yaw)
        yawMatrix[:, 1, 1] = np.cos(yaw)
        yawMatrix[:, 1, 2] = 0
        yawMatrix[:, 2, 0] = 0
        yawMatrix[:, 2, 1] = 0
        yawMatrix[:, 2, 2] = 1

        pitchMatrix[:, 0, 0] = np.cos(pitch)                                                                            # Pitch is rotation about the Y axis
        pitchMatrix[:, 0, 1] = 0
        pitchMatrix[:, 0, 2] = np.sin(pitch)
        pitchMatrix[:, 1, 0] = 0
        pitchMatrix[:, 1, 1] = 1
        pitchMatrix[:, 1, 2] = 0
        pitchMatrix[:, 2, 0] = -np.sin(pitch)
        pitchMatrix[:, 2, 1] = 0
        pitchMatrix[:, 2, 2] = np.cos(pitch)

        rollMatrix[:, 0, 0] = 1
        rollMatrix[:, 0, 1] = 0
        rollMatrix[:, 0, 2] = 0
        rollMatrix[:, 1, 0] = 0
        rollMatrix[:, 1, 1] = np.cos(roll)
        rollMatrix[:, 1, 2] = -np.sin(roll)
        rollMatrix[:, 2, 0] = 0
        rollMatrix[:, 2, 1] = np.sin(roll)
        rollMatrix[:, 2, 2] = np.cos(roll)

        if (not isarraybased):
            yawMatrix = yawMatrix.reshape([3, 3])
            pitchMatrix = pitchMatrix.reshape([3, 3])
            rollMatrix = rollMatrix.reshape([3, 3])

        R = yawMatrix @ pitchMatrix @ rollMatrix
        return RotationMatrix(R)

    # -----------------------------------------------------------------------------
    #           from_azimuth_elevation_roll
    # -----------------------------------------------------------------------------
    @staticmethod
    def from_azimuth_elevation_roll(azimuthradians: Union[np.ndarray, float],
                                    elevationradians: Union[np.ndarray, float],
                                    rollradians: Union[np.ndarray, float]) -> 'RotationMatrix':
        """
        Create a rotation matrix that rotates a vector by the specified azimuth, elevation and roll. It is is suitable
        for systems where the :math:`\\hat{x}` is horizontal and forwards, the :math:`\\hat{z}` axis is upwards and the
        :math:`\\hat{y}` points to the portside (left).  In this system, azimuth is a left handed (clockwise, compass) rotation
        around :math:`\\hat{z}`, elevation is a left handed rotation around :math:`\\hat{y}` (up is positive) and roll is a right-handed
        rotation around :math:`\\hat{x}` (leaning right is positive). The rotation matrix is created in the rotation order Azimuth,
        Elevation, Roll (`rzyx`)

        Parameters:
          azimuthradians
            The azimuthual rotation in radians. The rotation is a left-handed rotation around the :math:`\\hat{z}` axis.

          elevationradians
            The vertical rotation in radians. The rotation is a left-handed rotation around the :math:`\\hat{y}` axis.

          rollradians
            The roll rotation in radians. The rotation is in a right haded sense around the :math:`\\hat{x}` axis

        """

        return RotationMatrix.from_yaw_pitch_roll(-azimuthradians, -elevationradians, rollradians)

    # -----------------------------------------------------------------------------
    #           transform_vectors
    # -----------------------------------------------------------------------------

    def transform_vectors(self, lookvector: np.ndarray) -> np.ndarray:
        """
        Applies the internal rotation matrix to an array of N vectors.
        The incoming array must of the form (3,N) where the first dimension contains the X,Y,Z components of the N
        vectors.

        Parameters:
          lookvector:
            An array (3,N) of incoming vectors that will be rotated/transformed. The first dimension contains the X,Y,Z
            components of the N vectors

        Returns:
          An array of size (3,N) which is the set of vectors generated by applying the rotation matrix to the incoming vectors.
          The first dimension contains the X,Y,Z components of the rotated N vectors

        """

        v = self._R @ lookvector
        return v

    # -----------------------------------------------------------------------------
    #               get_yaw_pitch_roll
    # -----------------------------------------------------------------------------
    def yaw_pitch_roll(self) -> Tuple[float, float, float]:
        """
         Returns the current rotation matrix as the equivalent yaw, pitch roll (in radians).
         For example, see https://lavalle.pl/planning/node103.html

         Returns:
           A three element tuple of [yaw, pitch roll], all expressed in radians.
         """

        r = self.R
        r11 = r[0, 0]
        r21 = r[1, 0]
        r31 = r[2, 0]
        r32 = r[2, 1]
        r33 = r[2, 2]
        sinbeta = -r31  # get (3,1)
        cosbeta = np.sqrt(r32 * r32 + r33 * r33)  # r[3,2]^2 + r[3,3]^2
        pitch = np.arctan2(sinbeta, cosbeta)
        yaw = np.arctan2(r21, r11)
        roll = np.arctan2(r32, r33)
        return yaw, pitch, roll

    # -----------------------------------------------------------------------------
    #               from_ypr_in_first_frame_to_ypr_in_second_frame
    # -----------------------------------------------------------------------------
    @staticmethod
    def from_ypr_in_first_frame_to_ypr_in_second_frame(R1: 'RotationMatrix', RT: 'RotationMatrix') -> 'RotationMatrix':
        '''
        Converts from a yaw-pitch-roll, rotation matrix provided in one coordinate frame (Frame 1) to the equivalent rotation matrix in another coordinate frame (Frame 2).
        Returns the converted rotation matrix, :math:`R_2`, which should be used in the second coordinate frame. to apply the yaw-pitch-roll measured in the first frame

        This code was initially developed to convert and compare yaw, pitch roll measured on one IMU to yaw pitch and roll measured on another, rotated, IMU.
        It was applied to map the CNES IMU attitude solution to the iFTS IMU which was on the same gondola but at a rotated angle (nominally 90 degres). It was used to
        replace the iFTS IMU attitude solution which was at too coarse a resolution.

        **Theory**

        We use the rotation matrix defined for the first coordinate frame, :math:`R_1`, to map a vector, :math:`u_i`, specified in frame 1 coordinates
        to a geographic vector, :math:`w_{geo}`,

        ..  math::

                w_{geo} = R_1 u_i

        We wish to find the equivalent rotation matrix for frame 2, :math:`R_2` that maps the same vector but specified in frame 2, :math:`v_i`,
        to the same geographic vector, :math:`w_{geo}`. The vectors :math:`u_i` and :math:`v_i` are related via the rotation matrix transform from frame 2 to frame 1, :math:`R_T`.
        We wish to keep the vector given in frame 2, :math:`v_i`, constant and obtain its coordinates as expressed in frame 1. This is given using
        the inverse of the transform rotation matrix:

        ..  math::

            u_i = R_T^{-1} v_i

        thus we can substitute for :math:`u_i` in the above to get,

        ..  math::

            w_{geo} =  R_1 R_T^{-1} v_i

        Quite generally, the required rotation matrix in the second coordinate frame, :math:`R2` is the one that maps
        the instrument control frame to the platform vector,  it satisfies the equation,

        ..  math::

            w_{geo} = R_2 v_i

        Therefore, by inspection we can see the required rotation matrix in the second frame, :math:`R_2`, is given by:

        ..  math::

            R_2 = R_1 R_T^{-1}

        Parameters:
            R1:
                The rotation matrix in the first coordinate system, :math:`R_1`, For example, this may be derived from
                the yaw, pitch roll measured on by the first IMU system.

            RT:
                The rotation matrix that rotates from the second coordinate system to the first coordinate system, :math:`R_T`.

        Returns:
            RotationMatrix
                The rotation matrix, :math:`R2`, in the second coordinate frame.
        '''

        R2 = RotationMatrix()
        R2._R = R1.R @ RT.RInv
        return R2


# -----------------------------------------------------------------------------
#           class UnitVectors:
# -----------------------------------------------------------------------------
class UnitVectors:
    """
    A class to hold the unit vectors :math:`\\begin{bmatrix} \\hat{x} & \\hat{y} & \\hat{z} \\end{bmatrix}` for one
    reference frame specified in terms of another reference frame, :math:`\\begin{bmatrix} \\hat{x}' & \\hat{y}' & \\hat{z}' \\end{bmatrix}`
    The components of each unit vector occupy a column of the matrix.

        .. math::

             \\begin{pmatrix}
             \\hat{x}.\\hat{x}' & \\hat{y}.\\hat{x}' & \\hat{z}.\\hat{x}' \\\\
             \\hat{x}.\\hat{y}' & \\hat{y}.\\hat{y}' & \\hat{z}.\\hat{y}' \\\\
             \\hat{x}.\\hat{z}' & \\hat{y}.\\hat{z}' & \\hat{z}.\\hat{z}' \\\\
             \\end{pmatrix}

    Note that the components of :math:`\\begin{bmatrix} \\hat{x}' & \\hat{y}' & \\hat{z}' \\end{bmatrix}` in terms of
    :math:`\\begin{bmatrix} \\hat{x} & \\hat{y} & \\hat{z} \\end{bmatrix}` is the transpose (which is also the inverse)
    of this matrix.

    """
    def __init__(self, vectors: Union[List[List[float]], Tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray] = None):
        self._R: np.ndarray = RotationMatrix.IUnit() if vectors is None else np.array(vectors).transpose()

    @property
    def R(self):
        return self._R

    def X(self) -> np.ndarray:
        return self._R[:, 0]

    def Y(self) -> np.ndarray:
        return self._R[:, 1]

    def Z(self) -> np.ndarray:
        return self._R[:, 2]
