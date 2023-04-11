from abc import abstractmethod, ABC

class NeutrinoOscillation(ABC):
    """Abstract class to calculate neutrino oscillation"""

    def __init__(SNspectrum:sn_spectrum):
        self.sn_spectrum = sn_spectrum

    @abstractmethod
    def get_flux_nue():
        pass

    @abstractmethod
    def get_flux_anue():
        pass

    @abstractmethod
    def get_flux_nux():
        pass

    @abstractmethod
    def get_fluxes():
        pass


    
