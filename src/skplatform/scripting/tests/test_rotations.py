import numpy as np

def test_calculate_2d_rotation_matrices():
    '''
    Tests the correctness of the matrices used for rotating the
    central line-of-sight over a 2D field of view.
    '''
    from skplatform.scripting.rotations import calculate_2d_rotation_matrices

    results = calculate_2d_rotation_matrices(10.0, 6.0, 5, 5)
    # results should contain 25 rotation matrices
    # with azimuth angles -3, -1.5, 0, 1.5, 3
    # and elevation angles -5, -2.5, 0, 2.5, 5
    expected = np.array([
        [[0.99482945,-0.05233596,0.08703630],  # -3, -5
        [0.99585333,-0.02617695,0.08712588],  # -1.5, -5
        [0.99619470,0.00000000,0.08715574],  #  0, -5
        [0.99585333,0.02617695,0.08712588],  # 1.5, -5
        [0.99482945,0.05233596,0.08703630]],  #  3, -5

        [[0.99767906,-0.05233596,0.04355961],  # -3, -2.5
        [0.99870587,-0.02617695,0.04360444],  # -1.5, -2.5
        [0.99904822,0.00000000,0.04361939],  #  0, -2.5
        [0.99870587,0.02617695,0.04360444],  # 1.5, -2.5
        [0.99767906,0.05233596,0.04355961]],  #  3, -2.5

        [[0.99862953,-0.05233596,0.00000000],  # -3, -0
        [0.99965732,-0.02617695,0.00000000],  # -1.5, -0
        [1.00000000,0.00000000,0.00000000],  #  0, -0
        [0.99965732,0.02617695,0.00000000],  # 1.5, -0
        [0.99862953,0.05233596,0.00000000]],  #  3, -0

        [[0.99767906,-0.05233596,-0.04355961],  # -3,  2.5
        [0.99870587,-0.02617695,-0.04360444],  # -1.5,  2.5
        [0.99904822,0.00000000,-0.04361939],  #  0,  2.5
        [0.99870587,0.02617695,-0.04360444],  # 1.5,  2.5
        [0.99767906,0.05233596,-0.04355961]],  #  3,  2.5

        [[0.99482945,-0.05233596,-0.08703630],  # -3,  5
        [0.99585333,-0.02617695,-0.08712588],  # -1.5,  5
        [0.99619470,0.00000000,-0.08715574],  #  0,  5
        [0.99585333,0.02617695,-0.08712588],  # 1.5,  5
        [0.99482945,0.05233596,-0.08703630]],  #  3,  5
    ])
    x = np.array([1.0, 0.0, 0.0])
    for u in range(5):
        for v in range(5):
            assert np.allclose(results[u,v]@x, expected[u,v]), \
            f'Allclose X rotation failed: {u},{v}\t\n{results[u,v]@x}\t\n{expected[u,v]}'

    # test remaining vectors with fewer points
    results = calculate_2d_rotation_matrices(10.0, 6.0, 3, 3)

    # Y Axis
    expected = np.array([
        [[0.05213680,0.99862953,0.00456138],  # -3, -5
        [0.00000000,1.00000000,0.00000000],  # -0, -5
        [-0.05213680,0.99862953,-0.00456138]],  #  3, -5

        [[0.05233596,0.99862953,0.00000000],  # -3,  0
        [0.00000000,1.00000000,0.00000000],  # -0,  0
        [-0.05233596,0.99862953,0.00000000]],  #  3,  0

        [[0.05213680,0.99862953,-0.00456138],  # -3,  5
        [0.00000000,1.00000000,0.00000000],  # -0,  5
        [-0.05213680,0.99862953,0.00456138]],  #  3,  5
    ])

    y = np.array([0.0, 1.0, 0.0])
    for u in range(3):
        for v in range(3):
            assert np.allclose(results[u,v]@y, expected[u,v]), \
            f'Allclose Y rotation failed: {u},{v}\t\n{results[u,v]@y}\t\n{expected[u,v]}'

    # Z Axis
    expected = np.array([
        [[-0.08715574,0.00000000,0.99619470],  # -3, -5
        [-0.08715574,0.00000000,0.99619470],  # -0, -5
        [-0.08715574,0.00000000,0.99619470]],  #  3, -5

        [[0.00000000,0.00000000,1.00000000],  # -3,  0
        [0.00000000,0.00000000,1.00000000],  # -0,  0
        [0.00000000,0.00000000,1.00000000]],  #  3,  0

        [[0.08715574,0.00000000,0.99619470],  # -3,  5
        [0.08715574,0.00000000,0.99619470],  # -0,  5
        [0.08715574,0.00000000,0.99619470]],  #  3,  5
    ])
    
    z = np.array([0.0, 0.0, 1.0])
    for u in range(3):
        for v in range(3):
            assert np.allclose(results[u,v]@z, expected[u,v]), \
            f'Allclose Z rotation failed: {u},{v}\t\n{results[u,v]@z}\t\n{expected[u,v]}'

    # Arbitrary Vector
    expected = np.array([
        [[0.95981051,0.94629358,1.08779238],  # -3, -5
        [0.90903896,1.00000000,1.08335044],  # -0, -5
        [0.85553690,1.05096549,1.07866962]],  #  3, -5
        
        [[1.05096549,0.94629358,1.00000000],  # -3,  0
        [1.00000000,1.00000000,1.00000000],  # -0,  0
        [0.94629358,1.05096549,1.00000000]],  #  3,  0
        
        [[1.13412199,0.94629358,0.90459702],  # -3,  5
        [1.08335044,1.00000000,0.90903896],  # -0,  5
        [1.02984839,1.05096549,0.91371978]],  #  3,  5
    ])

    v = np.array([1.0, 1.0, 1.0])
    for i in range(3):
        for j in range(3):
            assert np.allclose(results[i,j]@v, expected[i,j]), \
            f'Allclose V rotation failed: {i},{j}\t\n{results[i,j]@v}\t\n{expected[i,j]}'
    
    # logging doesn't work unless use warning
    print(f'test_calculate_2d_rotation_matrices() passed.\n'
          f'\tRotation matrices for a 2d grid of pixels are correctly constructed.')


def plot_calculate_2d_rotation_matrices():
    import matplotlib
    matplotlib.use('Agg')
    from skplatform.scripting.visualization.plotlib3d import basic_setup, plot_vector
    from skplatform.scripting.rotations import calculate_2d_rotation_matrices, calculate_elevation_rotation_matrix

    R = calculate_2d_rotation_matrices(30.0, 20.0, 10, 10)
    fig, ax = basic_setup()

    # rvec = np.zeros((R.shape[0], R.shape[1], 3))
    x = np.array([1.0, 0.0, 0.0])
    for k in range(R.shape[0]):
        for j in range(R.shape[1]):
            rvec = R[k,j] @ x
            plot_vector(ax[0], rvec[0], rvec[1], rvec[2])
    ax[0].view_init(azim=45, elev=45)
    fig.savefig('tangent_2d_los.png')

    # nadir
    fig, ax = basic_setup()

    xR = calculate_elevation_rotation_matrix(90.0)
    R = calculate_2d_rotation_matrices(60.0, 80.0, 100, 100)
    # rvec = np.zeros((R.shape[0], R.shape[1], 3))
    x = np.array([1.0, 0.0, 0.0])
    x_len = 1.0 - 0.01
    for k in range(R.shape[0]):
        for j in range(R.shape[1]):
            rvec =  xR @ R[k,j] @ x
            plot_vector(ax[0], rvec[0], rvec[1], rvec[2],
                        rvec[0]*x_len, rvec[1]*x_len, rvec[2]*x_len,
                        )
    ax[0].view_init(azim=45, elev=45)
    fig.savefig('nadir_2d_los.png')

def test_construct_nadir_rotation_matrix():
    '''
    Tests the correctness of the matrix constructed to rotate 
    the lines-of-sight into a nadir orientation. This test relies on
    a correct calculate_2d_rotation_matrices() function.
    '''
    from scipy.spatial.transform import Rotation
    
    from skplatform.scripting.propagator import NadirPropagator
    from skplatform.scripting.rotations import calculate_2d_rotation_matrices

    fwd = np.array([1.0, 0.0, 0.0])
    rev = np.array([-1.0, 0.0, 0.0])
    left = np.array([0.0, 1.0, 0.0])
    right = np.array([0.0, -1.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    dwn = np.array([0.0, 0.0, -1.0])


    local_up = np.array([0.0, 0.0, 1.0])
    bore_sight = np.array([1.0, 0.0, 0.0])
    mat = NadirPropagator.construct_nadir_rotation_matrix(local_up, bore_sight)

    # should rotate a forward vector down, a reverse up
    assert np.allclose(mat@fwd, dwn), \
            f'{mat@fwd} {dwn}'
    assert np.allclose(mat@rev, up)
    # shouldn't rotate a sideways vector 
    assert np.allclose(mat@left, left)
    assert np.allclose(mat@right, right)
    # z axis should become x
    assert np.allclose(mat@up, fwd)
    assert np.allclose(mat@dwn, rev)

    
    local_up = np.array([-0.28390867, -0.24403605,  0.9272768 ])
    bore_sight = np.array([-0.88058006, -0.31633352, -0.35286239])
    sideways = np.cross(bore_sight, -local_up)
    # nadir matrix maps forward to down, across to across, and up to along track
    N = NadirPropagator.construct_nadir_rotation_matrix(local_up, bore_sight)

    R = calculate_2d_rotation_matrices(180.0, 0, 3, 1)  # forward, along track, elevation
    # fwd ->R-> up ->N-> fwd
    assert np.allclose(N@R[0,0]@fwd, bore_sight), f'{N@R[0,0]@fwd}, {bore_sight}' # back 90 around y, then fwd 90
    assert np.allclose(N@R[1,0]@fwd, -local_up), f'{N@R[0,2]@fwd}, {-local_up}'
    assert np.allclose(N@R[2,0]@fwd, -bore_sight), f'{N@R[0,2]@fwd}, {-bore_sight}' # fwd 90 around y, then another 90
    
    R = calculate_2d_rotation_matrices(0, 180.0, 1, 3)  # sideways, across track, azimuth
    assert np.allclose(N@R[0,0]@fwd, -sideways), f'{N@R[0,0]@fwd}, {-sideways}'
    #assert np.allclose(N@R[0,0]@fwd, -sideways), f'{N@R[0,0]@fwd}, {-sideways}'
    assert np.allclose(N@R[0,2]@fwd, sideways), f'{N@R[0,2]@fwd}, {sideways}'

    R = calculate_2d_rotation_matrices(0, 90.0, 1, 3)  # +/- 45 sideways, across track, azimuth
    r2 = np.sin(np.pi/4)  # == cos()
    # this is rotated 45 to the right by R so becomes root2*x, -root2*y
    # then N rotated the frame so +x becomes -z, y unchanged -> -root2*y, -root2*z
    # y is sideways, z is down (-local_up)
    assert np.allclose(N@R[0,0]@fwd, -sideways*r2 + -local_up*r2), f'{N@R[0,0]@fwd}, {-sideways*r2 + -local_up*r2}'
    assert np.allclose(N@R[0,2]@fwd, sideways*r2 + -local_up*r2), f'{N@R[0,2]@fwd}, {sideways*r2 + -local_up*r2}'

    R = calculate_2d_rotation_matrices(90.0, 0.0, 3, 1)  # +/- 45 sideways, across track, azimuth
    # this is rotated 45 to the right by R so becomes root2*x, -root2*y
    # then N rotated the frame so +x becomes -z, y unchanged -> -root2*y, -root2*z
    # y is sideways, z is down (-local_up)
    assert np.allclose(N@R[0,0]@fwd, bore_sight*r2 + -local_up*r2), f'{N@R[0,0]@fwd}, {bore_sight*r2 + -local_up*r2}'
    assert np.allclose(N@R[2,0]@fwd, -bore_sight*r2 + -local_up*r2), f'{N@R[2,0]@fwd}, {-bore_sight*r2 + -local_up*r2}'

    nR = Rotation.from_matrix(N)
    rR = Rotation.from_euler('Y', [[45.0]], degrees=True)
    # rR = Rotation.from_euler('zy', [[45.0, 45.0]], degrees=True)
    # rR = Rotation.from_euler('YZ', [[45.0, 45.0]], degrees=True)
    # rR = Rotation.from_euler('yz', [[45.0, 45.0]], degrees=True)
    assert np.allclose(N@R[2,0]@fwd, nR.as_matrix() @ rR.as_matrix() @ fwd),\
        f'{N @ R[2,0]@fwd} {nR.as_matrix() @ rR.as_matrix() @ fwd}'
    
    # reaching this point means all asserts passed
    print('test_construct_nadir_rotation_matrix() passed.\n'
          '\tThe matrix for rotating LOS to nadir is correct.')
    

if __name__ == "__main__":
    test_calculate_2d_rotation_matrices()
    test_construct_nadir_rotation_matrix()
    # #plot_calculate_2d_rotation_matrices()
    