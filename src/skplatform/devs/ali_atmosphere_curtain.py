import xarray as xr
import numpy as np
import pandas as pd
from skplatform.platform import OpticalGeometry
from skplatform.geodesy import Geodetic


class SimulationAtmosphere:

    def __init__(self, filename, angular_resoluiton=0.5, vertical_resolution=0.5):

        self._file = filename
        self.angular_resolution = angular_resoluiton
        self.vertical_resolution = vertical_resolution
        self.max_altitude = 45
        self.min_altitude = 0
        self.max_angle = 90
        self.min_angle = 0

        self.particle_size_2d = False
        self.cloud_effective_radius = 70

    @property
    def orbit_angle(self):
        return self.orbit_angle_bins[0:-1] + np.diff(self.orbit_angle_bins) / 2

    @property
    def altitude(self):
        return self.altitude_bins[0:-1] + np.diff(self.altitude_bins) / 2

    @property
    def orbit_angle_bins(self):
        return np.arange(self.min_angle, self.max_angle + self.angular_resolution / 2, self.angular_resolution)

    @property
    def altitude_bins(self):
        return np.arange(self.min_altitude, self.max_altitude + self.vertical_resolution / 2, self.vertical_resolution)

    @property
    def normal_vector(self):
        geo = Geodetic()
        l0 = geo.xyz_from_llh([self.latitude(0.0), self.longitude(0.0), 0.0])
        l1 = geo.xyz_from_llh([self.latitude(1.0), self.longitude(1.0), 0.0])
        orbit_dir = l1 - l0
        orbit_dir /= np.linalg.norm(orbit_dir)
        return np.cross(orbit_dir, l0 / np.linalg.norm(l0))

    @property
    def reference_vector(self):
        geo = Geodetic()
        location = geo.xyz_from_llh([self.latitude(0.0), self.longitude(0.0), 0.0])
        return location / np.linalg.norm(location)

    @property
    def aerosol(self):

        rg = np.ones((len(self.orbit_angle), len(self.altitude))) * 0.08
        sg = np.ones((len(self.orbit_angle), len(self.altitude))) * 1.6
        particle_size = {'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': sg,
                         'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': rg}
        ps_clim = sk.ClimatologyUserDefined2D(self.orbit_angle, self.altitude * 1000, particle_size,
                                              self.reference_vector, self.normal_vector)
        opt_prop = sk.MieAerosol(ps_clim, 'H2SO4')

        if self.particle_size_2d:
            xsec = np.ones((len(self.orbit_angle), len(self.altitude)))
            for angle_idx, angle in enumerate(self.orbit_angle):
                lat = self.latitude(angle)
                lon = self.longitude(angle)
                for idx, alt in enumerate(self.altitude):
                    xsec[angle_idx, idx] = opt_prop.calculate_cross_sections(sk.MSIS90(),
                                                                             latitude=lat, longitude=lon,
                                                                             altitude=alt * 1000, mjd=53000.0,
                                                                             wavelengths=750.0).total[0]
            xsec = xr.DataArray(xsec, coords=[self.orbit_angle, self.altitude], dims=['angle', 'altitude'])
        else:
            xsec = opt_prop.calculate_cross_sections(sk.MSIS90(),
                                                     latitude=self.latitude(10), longitude=self.longitude(10),
                                                     altitude=1000.0, mjd=53000.0,
                                                     wavelengths=750.0).total[0]
        extinction = self._downsample_omps()
        aerosol = extinction / xsec * 1e-5
        aerosol = aerosol.where(np.isfinite(aerosol)).fillna(0.0).transpose('angle', 'altitude')
        aer_clim = sk.ClimatologyUserDefined2D(self.orbit_angle, self.altitude * 1000,
                                               {'SKCLIMATOLOGY_AEROSOL_CM3': aerosol.values},
                                               self.reference_vector, self.normal_vector)
        return sk.Species(opt_prop, aer_clim, species='SKCLIMATOLOGY_AEROSOL_CM3')

    @property
    def cloud(self):

        cloud = self._downsample_calipso()
        opt_prop = sk.BaumIceCrystal(70)
        xsec = opt_prop.calculate_cross_sections(sk.MSIS90(), latitude=self.latitude(10), longitude=self.longitude(10),
                                                 altitude=1000.0, mjd=53000.0, wavelengths=750.0).total[0]
        p = opt_prop.calculate_phase_matrix(sk.MSIS90(), latitude=self.latitude(10), longitude=self.longitude(10),
                                            altitude=1000.0, mjd=53000.0, wavelengths=[750.0],
                                            cosine_scatter_angles=[-1])[0, 0]
        bs_to_ext = (1 / p[0, 0]) * np.pi * 4
        nd = cloud * bs_to_ext / xsec * 1e-5
        nd = nd.fillna(0.0)
        cloud_clim = sk.ClimatologyUserDefined2D(self.orbit_angle, self.altitude * 1000,
                                                 {'icecloud': nd.transpose('angle', 'altitude').values},
                                                 self.reference_vector, self.normal_vector)
        return sk.Species(opt_prop, cloud_clim)

    @property
    def water_vapour(self):
        h2o = self._downsample_h2o()
        opt_prop = sk.HITRANChemical('h2o')
        clim = sk.ClimatologyUserDefined2D(self.orbit_angle, self.altitude * 1000,
                                           {'h2o': h2o.transpose('angle', 'altitude').values},
                                           self.reference_vector, self.normal_vector)
        return sk.Species(opt_prop, clim)

    @property
    def air(self):

        return sk.Species(sk.Rayleigh(), sk.MSIS90())

    def optical_geometry(self, orbit_angle, altitude=10.0, satellite_altitude=500.0):

        lat = self.latitude(orbit_angle)
        lon = self.longitude(orbit_angle)
        time = self.time(orbit_angle)
        mjd = (time - pd.Timestamp('1858-11-17')) / pd.Timedelta(1, 'D')
        tp = Geodetic()
        location = tp.xyz_from_llh([lat, lon, altitude * 1000])
        normal = self.normal_vector

        look = np.cross(location / np.linalg.norm(location), normal)
        re = 6372
        tp_to_sat = np.sqrt((re + satellite_altitude) ** 2 - (re + altitude) ** 2) * 1000
        obs = location - look * tp_to_sat
        hor = np.cross(look, obs / np.linalg.norm(obs))
        up = np.cross(hor, look)

        return OpticalGeometry(observer=obs, look_vector=look, local_up=up, mjd=mjd)

    def sastkran_atmosphere(self, aerosol=True, cloud=True, h2o=True):

        sk_atmo = sk.Atmosphere()
        sk_atmo['air'] = self.air

        if aerosol:
            sk_atmo['aerosol'] = self.aerosol

        if cloud:
            sk_atmo['cloud'] = self.cloud

        if h2o:
            sk_atmo['h2o'] = self.water_vapour

        return sk_atmo

    def _downsample_h2o(self):
        era = xr.open_dataset(self._file, group='ERA5')
        angle_along_orbit = self._angle_along_orbit(np.asarray(era.latitude.values, dtype=float),
                                                    np.asarray(era.longitude.values, dtype=float))
        nd = (era.pressure * 100) / (era.temperature * 1.38064852e-23) / (100 ** 3)
        m_air = 28.9647
        m_h2o = 18.02
        h2o = era.specific_humidity * nd * (m_air / m_h2o)
        h2o = xr.DataArray(h2o.where(h2o > 0).transpose('altitude', 'time').values,
                           coords=(era.altitude.values, angle_along_orbit),
                           dims=['altitude', 'angle']).fillna(0.0)

        # TODO: extend h2o clim above 40km using MSIS90?

        h2o = self._downsample_to_atmo_grid(h2o).fillna(0.0)
        return h2o

    def latitude(self, orbit_angle):
        lat, lon, angles = self._calipso_position()
        return np.interp(orbit_angle, angles, lat)

    def longitude(self, orbit_angle):
        lat, lon, angles = self._calipso_position()
        return np.interp(orbit_angle, angles, lon)

    def time(self, orbit_angle):
        mjd = self.mjd(orbit_angle)
        return pd.Timedelta(mjd, 'D') + pd.Timestamp('1858-11-17')

    def mjd(self, orbit_angle):
        calipso = xr.open_dataset(self._file, group='CALIPSO')
        lat, lon, angles = self._calipso_position()
        time = calipso.time.values

        mjds = (time - np.datetime64('1858-11-17')) / np.timedelta64(1, 'D')
        mjd = np.interp(orbit_angle, angles, mjds)
        return mjd

    def _calipso_position(self):
        calipso = xr.open_dataset(self._file, group='CALIPSO')
        latitude = np.asarray(calipso.latitude.values, dtype=float)
        longitude = np.asarray(calipso.longitude.values, dtype=float)
        angles = self._angle_along_orbit(latitude, longitude, 0.0)
        return latitude, longitude, angles

    def _downsample_to_atmo_grid(self, array):

        array = array.groupby_bins(array.altitude, bins=self.altitude_bins).mean()
        array = array.groupby_bins(array.angle, bins=self.orbit_angle_bins).mean()
        array['altitude_bins'] = self.altitude
        array['angle_bins'] = self.orbit_angle
        array = array.sortby('angle_bins').sortby('altitude_bins') \
                     .rename({'altitude_bins': 'altitude', 'angle_bins': 'angle'})
        return array

    def _angle_along_orbit(self, latitude, longitude, zero_lat=0.0):

        lat0 = zero_lat
        lon0 = np.interp(zero_lat, latitude, longitude)
        geo = Geodetic()
        location = geo.xyz_from_llh([lat0, lon0, 0.0])
        xyz0 = location / np.linalg.norm(location)
        angle_along_orbit = np.ones_like(latitude)
        for idx, (lat, lon) in enumerate(zip(latitude, longitude)):
            location = geo.xyz_from_llh([lat, lon, 0.0])
            xyz = location / np.linalg.norm(location)
            if lat < lat0:
                angle_along_orbit[idx] = -np.arccos(np.dot(xyz, xyz0))
            else:
                angle_along_orbit[idx] = np.arccos(np.dot(xyz, xyz0))

        return angle_along_orbit * 180 / np.pi

    def _downsample_calipso(self):

        calipso = xr.open_dataset(self._file, group='CALIPSO')
        angle_along_orbit = self._angle_along_orbit(np.asarray(calipso.latitude.values, dtype=float),
                                                    np.asarray(calipso.longitude.values, dtype=float))
        cloud = xr.DataArray(calipso.cloud_backscatter.where(calipso.cloud_backscatter > 0).values,
                             coords=(calipso.altitude.values, angle_along_orbit),
                             dims=['altitude', 'angle']).fillna(0.0)
        calipso_cloud = self._downsample_to_atmo_grid(cloud)

        return calipso_cloud

    def _downsample_omps(self):
        omps = xr.open_dataset(self._file, group='OMPS')
        angle_along_orbit = self._angle_along_orbit(np.asarray(omps.latitude.values, dtype=float),
                                                    np.asarray(omps.longitude.values, dtype=float))

        aerosol = xr.DataArray(omps.extinction.values.T,
                               coords=(omps.altitude.values, angle_along_orbit),
                               dims=['altitude', 'angle'])
        tropopause = xr.DataArray(omps.tropopause_altitude.values.T,
                                  coords=[angle_along_orbit],
                                  dims=['angle'])

        aerosol = aerosol.interp(altitude=self.altitude).interp(angle=self.orbit_angle)
        tropopause = tropopause.interp(angle=self.orbit_angle)

        # fill nans then decay aerosol below tropopause and above upper bound
        min_altitude = ((aerosol > 0) * aerosol.altitude)
        min_altitude = min_altitude.where(min_altitude > 0).min(dim='altitude').values
        max_altitude = ((aerosol > 0) * aerosol.altitude)
        max_altitude = max_altitude.where(max_altitude > 0).max(dim='altitude').values
        alts = aerosol.altitude.values
        for profile, min_alt, max_alt, trop in zip(aerosol.transpose('angle', 'altitude').values,
                                                   min_altitude, max_altitude, tropopause.values):
            min_idx = np.where(alts >= min_alt)[0][0]
            max_idx = np.where(alts <= max_alt)[0][-1]
            profile[alts < min_alt] = profile[min_idx]
            scale = np.exp(-1 * alts)
            scale *= (profile[max_idx] / scale[max_idx])
            profile[alts > max_alt] = scale[alts > max_alt]
            alt_below_trop = alts - trop
            alt_below_trop[alt_below_trop > 0] = 0
            alt_below_trop[alt_below_trop < -3] = -3
            alt_below_trop *= (1 / 3)
            alt_below_trop += 1
            profile *= alt_below_trop

        return aerosol
