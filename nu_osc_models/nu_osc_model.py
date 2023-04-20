from abc import abstractmethod, ABC
from .. import sn_spectra

class nu_osc_model(ABC):
    """Abstract class for calculating neutrino oscillation"""

    def __init__(self, SNspectra:sn_spectra):
        pass

    @abstractmethod
    def get_sn_spectra(self):
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
        pass

    @abstractmethod
    def get_number_spectra_anue(self):
        pass

    @abstractmethod
    def get_number_spectra_nux(self):
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
