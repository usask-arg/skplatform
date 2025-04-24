import numpy as np
# from numba import njit


# @njit(fastmath=True, cache=True)
def calculate_azimuth_rotation_matrix(angle):
    R = np.zeros((3, 3))
    rot_cos = np.cos(np.deg2rad(angle))
    rot_sin = np.sin(np.deg2rad(angle))
    R[0, 0] = rot_cos
    R[1, 1] = rot_cos
    R[0, 1] = -rot_sin
    R[1, 0] = rot_sin
    R[2, 2] = 1.0
    return R  # -> [xyz, xyz']


# @njit(fastmath=True, cache=True)
def calculate_elevation_rotation_matrix(angle):
    R = np.zeros((3, 3))
    rot_cos = np.cos(np.deg2rad(angle))
    rot_sin = np.sin(np.deg2rad(angle))
    R[0, 0] = rot_cos
    R[2, 0] = -rot_sin
    R[1, 1] = 1.0
    R[0, 2] = rot_sin
    R[2, 2] = rot_cos
    return R  # -> [xyz, xyz']


# @njit(fastmath=True, cache=True)
def calculate_roll_rotation_matrix(angle):
    R = np.zeros((3, 3))
    rot_cos = np.cos(np.deg2rad(angle))
    rot_sin = np.sin(np.deg2rad(angle))
    R[0, 0] = 1.0
    R[1, 1] = rot_cos
    R[1, 2] = -rot_sin
    R[2, 1] = rot_sin
    R[2, 2] = rot_cos
    return R  # -> [xyz, xyz']


# @njit(fastmath=True, cache=True)
def calculate_2d_rotation_matrices(along_fov: float, across_fov: float, along_pixels: int, across_pixels: int):
    """
    Calculates the rotation matrices for rotating the boresight to each of
    the pixels in the detector. Rotations are in ZYX axis order, or by
    azimuth angle and then by elevation angle. The matrices do not
    include a roll rotation.

    Angles are the full field-of-view. The rotations span from -fov/2 to
    +fov/2 with the end points included.

    Parameters
    ----------
    along_fov : float
        The full angle of the vertical field of view, the elevation angle.
        By default this is aligned with the satellite ground track in
        both limb and nadir geometries.
    across_fov : float
        The full angle of the horizontal field of view, the azimuth angle.
        By default this is perpendicular to the satellite ground track in
        both limb and nadir geometries.
    along_pixels : int
        The number of pixels in the vertical FOV.
    across_pixels : int
        The number of pixels in the horizontal FOV.
    """

    elv_angles = np.deg2rad(np.linspace(-along_fov/2, along_fov/2, along_pixels))
    azm_angles = np.deg2rad(np.linspace(-across_fov/2, across_fov/2, across_pixels))
    R = np.zeros((along_pixels, across_pixels, 3, 3))

    # elevation, rotation about the x axis, look vector looks futher or closer along boresight
    for elv, elevation_angle in enumerate(elv_angles):
        elv_R = np.zeros((3,3))
        elv_cos = np.cos(elevation_angle)
        elv_sin = np.sin(elevation_angle)
        elv_R[0, 0] =  elv_cos
        elv_R[2, 0] = -elv_sin
        elv_R[1, 1] = 1.0
        elv_R[0, 2] =  elv_sin
        elv_R[2, 2] =  elv_cos

        # azimuth, rotation about the z axis
        for azm, azimuth_angle in enumerate(azm_angles):
            azm_R = np.zeros((3,3))
            azm_cos = np.cos(azimuth_angle)
            azm_sin = np.sin(azimuth_angle)

            azm_R[0, 0] = azm_cos
            azm_R[1, 1] = azm_cos
            azm_R[0, 1] = -azm_sin
            azm_R[1, 0] = azm_sin
            azm_R[2, 2] = 1.0

            # print(f'{np.rad2deg(xa)} {np.rad2deg(ya)}')
            R[elv, azm, :, :] = elv_R @ azm_R
    return R
