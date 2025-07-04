import unittest
import numpy as np
from skplatform import Platform


class SkretrievalPlatformTests2(unittest.TestCase):

    def test_instrument_mounting(self):
        """
        Tests the instrument mounting in the platform control frame
        """
        platform = Platform()
        vicf = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1))).transpose()
        lat = np.linspace(0.0, 90.0, 1000)
        lng = np.linspace(0.0, 180.0, 1000)
        height = np.full([1000], 650000.0)
        position = np.array((lat, lng, height)).transpose()
        platform.acs.set_platform_location(latlonheight=position)

        omegax = np.radians(np.linspace(-2.5, 2.5, 500))
        omegay = np.radians(np.linspace(-1.2, 0.0, 400))
        omega = np.zeros([omegay.size, omegax.size, 2])
        for ix in range(omegax.size):
            omega[:, ix, 0] = omegay
        for iy in range(omegay.size):
            omega[iy, :, 1] = omegax

        cosomegay = np.cos(omega[..., 0])
        vicf = np.zeros( list(omega.shape[:-1]) + [3])
        vicf[..., 0] = cosomegay * np.cos(omega[..., 1])
        vicf[..., 1] = cosomegay * np.sin(omega[..., 1])
        vicf[..., 2] = np.sin(omega[..., 0])

        platform.acs.mount_instrument_on_platform(60, 30.0, 60.0)
        vgeo = platform.acs.convert_icf_to_gcf(vicf)
        vecef = platform.acs.convert_icf_to_ecef(vicf)
        vecefsmall = platform.acs.convert_icf_to_ecef(vicf[10,20,:])

        platform.acs.set_platform_location(latlonheight=position[50,:])
        vgeo = platform.acs.convert_icf_to_gcf(vicf)
        vecef = platform.acs.convert_icf_to_ecef(vicf)
        vecefsmall = platform.acs.convert_icf_to_ecef(vicf[10,20,:])

        # print('rotation = (60,30,0)')
        # print('Geographic Unit vectors')
        # print(vgeo)
        # print('ECEF Unit vectors')
        # print(vecef)



if __name__ == "__main__":
    tests = SkretrievalPlatformTests2()
    unittest.main()
