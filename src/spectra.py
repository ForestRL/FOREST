from abc import ABCMeta, abstractmethod, ABC


class spectra(metaclass=ABCMeta):
    """ addresses neutrino spectra
    """

    @abstractmethod
    def get_luminosities(self):
        """calculates luminosities for each flavor.
            return:
                numpy arrays
                time, nue lum, anue lum, nux_lum
        """
        pass

    @abstractmethod
    def get_fluxes(self):
        """calculates fluxes for each flavor.
            return:
                numpy arrays
                time, nue flux, anue flux, nux_flux
        """
        pass 

    @abstractmethod
    def get_average_energies(self):
        """ calculates average energies for each flavor.
            the formula in the Nakazato format is
            Luminsity/Flux/1.6022e-6
            return:
                numpy arrays
                time, nue ene, anue ene, nux ene
        """
        pass

    @abstractmethod
    def get_ene1_bins(self):
        """ returns the array of ene1_bins.
            return:
                   numpy array
                   ene1_bins
        """
        pass

    @abstractmethod
    def get_ene2_bins(self):
        """ returns the array of ene2_bins.
            return:
                   numpy array
                   ene2_bins
        """
        pass

    @abstractmethod
    def get_ene_center_bins(self):
        """ returns the array of ene_center_bins.
            return:
                   numpy array
                   ene_center_bins
        """
        pass
    
    @abstractmethod
    def get_nue_energy_spectra(self):
        """ returns the array of nue_energy_spectra.
            return:
                   numpy array
                   
        """
        pass


    @abstractmethod
    def get_anue_energy_spectra(self):
        """ returns the array of anue_energy_spectra.
            return:
                   numpy array
                   
        """
        pass

    @abstractmethod    
    def get_nux_energy_spectra(self):
        """ returns the array of nux_energy_spectra.
            return:
                   numpy array
                   
        """
        pass

    @abstractmethod
    def get_nue_number_spectra(self):
        """ returns the array of nue_number_spectra.
            return:
                   numpy array
                   
        """
        pass
    
    @abstractmethod
    def get_anue_number_spectra(self):
        """ returns the array of anue_number_spectra.
            return:
                   numpy array
                   
        """
        pass

    @abstractmethod
    def get_nux_number_spectra(self):
        """ returns the array of nux_number_spectra.
            return:
                   numpy array
                   
        """
        pass

    @abstractmethod
    def get_times(self):
        """ returns the array of times.
            return:
                   nunpy array
        """
        pass

    @abstractmethod
    def get_distance(self) -> float:
        """ returns the distance of the SN from earth.
            return:
                   distance (kpc)
        """
        pass
