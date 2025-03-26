import numpy as np
import scipy.signal


# -----------------------------------------------------------------------------
#               SolarSpectrum
# -----------------------------------------------------------------------------
class SolarSpectrum:
    '''
    A class used to load solar irradiance spectrum. The class currently supports the following spectra:

    #. SAO2010
    #. Fontela UVIS to 3 microns
    #. Fontela UVIS to 100 micron,

    The choice of solar spectrum is specified in the constructor, the default is `sao2010`. The other values
    are 'fontela_uvis3micron` and `fontela_uvis100micron'.

    The solar data are implemented as python data arrays stored in the source code and can a second or two to load. But
    we have written the code so the solar data are only imported in for the model you request and only when you request it.
    '''

    def __init__(self, spectrum_name: str = 'sao2010', fwhm=None):
        self._solardistanceinAU = 1.0
        self._samplespacing = 0.01
        self._resolution_nm = 0.05
        self._wavelen_nm: np.ndarray = None
        self._irradiance: np.ndarray = None
        self._loadspectrum(spectrum_name)
        if fwhm is not None:
            smooth = self._smooth_spectrum(fwhm)
            self._irradiance = smooth
            self._resolution_nm = fwhm

    def _load_sao2010(self):
        from .solar_spectrum_impl.sao2010 import solar_spectrum_sao2010
        self._wavelen_nm, self._irradiance = solar_spectrum_sao2010()
        self._samplespacing = 0.01
        self._resolution_nm = 0.04  # FWHM in nm

    def _load_fontela_3micron(self):
        from .solar_spectrum_impl.fontela_uvis3micron import solar_spectrum_fontelaspectrum
        self._wavelen_nm, self._irradiance = solar_spectrum_fontelaspectrum()
        self._samplespacing = 0.02
        self._resolution_nm = 0.1  # FWHM in nm

    def _load_fontela_100micron(self):
        from .solar_spectrum_impl.fontela_uvis100micron import solar_spectrum_fontelaspectrum100
        self._wavelen_nm, self._irradiance = solar_spectrum_fontelaspectrum100()
        self._samplespacing = 0.2
        self._resolution_nm = 1.0           # FWHM in nm

    def _loadspectrum(self, spectrum_name: str):

        name = spectrum_name.lower()
        if (name == 'sao2010'):
            self._load_sao2010()
        elif (name == 'fontela_uvis3micron'):
            self._load_fontela_3micron()
        elif (name == 'fontela_uvis100micron'):
            self._load_fontela_100micron()
        else:
            raise ValueError("Invalid solar spectrum value {}. Valid values for the solar spectrum are 'sao2010', 'fontela_uvis3micron', or 'fontela_uvis100micron'".format(spectrum_name))

        # -----------------------------------------------------------------------------
        #               irradiance
        # -----------------------------------------------------------------------------
    def irradiance(self, nm: np.ndarray):
        '''
        Returns the solar irradiance at 1 A.U. at the requested wavelengths.

        Parameters
        ----------
        nm : np.ndarray
            An array of wavelengths in nanometers. The current solar irradiance table is linearly interpolated to
            the desired wavelengths. The irradiance signal is set to 0.0 for all wavelengths outside of the current
            solar table wavelength bounds.

        Returns
        -------
        np.ndarray
            The solar irradiance for each of the requested wavelengtths
        '''
        irrad = np.interp(nm, self._wavelen_nm, self._irradiance, left=0.0, right=0.0) / (self._solardistanceinAU * self._solardistanceinAU)
        return irrad

    # -----------------------------------------------------------------------------
    #               set_solar_distance_au
    # -----------------------------------------------------------------------------
    def set_solar_distance_au(self, au: float):
        '''
        Sets the distance of the sun in astronomical units (AU).

        Parameters
        ----------
        au: float
            The distance of the Sun in astronomical units. The default upon initialization is 1.0

        '''
        self._solardistanceinAU = au

    @property
    def sample_spacing(self) -> float:
        '''
        The spacing in nanometers between spectral samples in the table.
        '''
        return self._samplespacing

    @property
    def resolution_fwhm(self) -> float:
        '''
        The FWHM in nanometers of the instrument resolution used to measure the current solar spectrum.
        '''
        return self._resolution_nm

    @property
    def min_valid_wavelength(self):
        '''
        The minimum wavelength in nanometers in the current table.
        '''
        return self._wavelen_nm[0]

    @property
    def max_valid_wavelength(self):
        '''
        The maximum wavelength in nanometers in the curreent table
        '''
        return self._wavelen_nm[-1]

    def _smooth_spectrum(self, fwhm_nm) -> np.ndarray:

        stdev = fwhm_nm / 2.3548200450309493820231386529194                             # Get the stanradr deviation required by the user
        hires_stdev = self.resolution_fwhm / 2.3548200450309493820231386529194          # get the standard deviation of the high resolution solar signal
        sigma2 = (stdev ** 2 - hires_stdev ** 2)                                        # get the modified standrad deviation, subtract the solar standard deviation
        assert sigma2 > 0.0
        sigma = np.sqrt(sigma2)
        dw = self.sample_spacing
        x = np.arange(-20 * sigma, 20 * sigma, dw)
        g = np.exp(-0.5 * np.square(x) / sigma2) / (np.sqrt(2.0 * np.pi * sigma2)) * dw
        smooth = scipy.signal.fftconvolve(self._irradiance, g, mode='same')
        return smooth
