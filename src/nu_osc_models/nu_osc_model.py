import sys
sys.path.append("../")
from abc import ABCMeta, abstractmethod, ABC
from sn_spectra import SNspectra

class nu_osc_model(metaclass=ABCMeta):
    """Abstract class for calculating neutrino oscillation"""

    def __init__(self, sn_spectra:SNspectra):
        pass

    @abstractmethod
    def get_sn_spectra(self) -> SNspectra:
        """Returns SNspectra class"""
        pass

    @abstractmethod
    def get_flux_nue(self):
        pass

    @abstractmethod
    def get_flux_anue(self):
        pass

    @abstractmethod
    def get_flux_nux(self):
        pass

    @abstractmethod
    def get_fluxes(self):
        pass
    
    @abstractmethod
    def get_fluxes(self):
        pass

    @abstractmethod
    def get_lum_nue(self):
        pass
    
    @abstractmethod
    def get_lum_anue(self):
        pass

    @abstractmethod
    def get_lum_nux(self):
        pass

    @abstractmethod
    def get_lums(self):
        pass

    @abstractmethod
    def get_aveene_nue(self):
        pass
    
    @abstractmethod
    def get_aveene_anue(self):
        pass

    @abstractmethod
    def get_aveene_nux(self):
        pass

    @abstractmethod
    def get_aveenes(self):
        pass

    @abstractmethod
    def get_number_spectra_nue(self):
        """ returns the array of nue_number_spectra.
            return:
                   numpy array
        """
        pass

    @abstractmethod
    def get_number_spectra_anue(self):
        """ returns the array of anue_number_spectra.
            return:
                   numpy array                   
        """
        pass

    @abstractmethod
    def get_number_spectra_nux(self):
        """ returns the array of nux_number_spectra.
            return:
                   numpy array        
        """
        pass

    @abstractmethod
    def get_energy_spectra_nue(self):
        pass

    @abstractmethod
    def get_energy_spectra_anue(self):
        pass

    @abstractmethod
    def get_energy_spectra_nux(self):
        pass
