import unittest
import numpy as np
from skplatform import RotationMatrix


class MatrixTests(unittest.TestCase):

    def test_matrices(self):

        M = RotationMatrix.from_yaw_pitch_roll(np.radians(45), np.radians(30), np.radians(45))
        x = np.array([np.cos(np.radians(30)), 0.0, np.sin(np.radians(30))])
        y = M.R @ x
