import sys
sys.path.append("../")
from nu_osc_model import nu_osc_model
from sn_spectra import SNspectra

class MSW_normal(nu_osc_model):
    """Neutrino oscillationi class of the MSW effect with the normal hierarchy"""

    def __init__(self, sn_spectra:SNspectra):

        p = 1.0#0.69

        self.__sn_spectra = sn_spectra
        self.__fluxes_nue   = sn_spectra.get_nue_number_spectra() # array of flux of nue
        self.__fluxes_anue  = sn_spectra.get_anue_number_spectra() # array of flux of anue
        self.__fluxes_nux   = sn_spectra.get_nux_number_spectra() # array of flux of nux
        self.__lums_nue   = sn_spectra.get_nue_energy_spectra() # array of energy flux of nue
        self.__lums_anue  = sn_spectra.get_anue_energy_spectra() # array of energy flux of anue
        self.__lums_nux   = sn_spectra.get_nux_energy_spectra() # array of energy flux of nux


        # neutrino number spectra after oscillation
        self.__osc_fluxes_nue = self.__fluxes_nux 
        self.__osc_fluxes_anue = p*self.__fluxes_anue + (1-p)*self.__fluxes_nux
        self.__osc_fluxes_nux = (self.__fluxes_nue + (1-p)*self.__fluxes_anue + (2+p)*self.__fluxes_nux)*0.25

        # neutrino energy spectra after oscillation
        self.__osc_lums_nue = self.__lums_nux
        self.__osc_lums_anue = p*self.__lums_anue + (1-p)*self.__lums_nux
        self.__osc_lums_nux = (self.__lums_nue + (1-p)*self.__lums_anue + (2+p)*self.__lums_nux)*0.25


    def get_sn_spectra(self):
        return self.__sn_spectra

    
    def get_flux_nue():
        pass

    
    def get_flux_anue():
        pass

    
    def get_flux_nux():
        pass


    def get_fluxes():
        pass
    
    
    def get_fluxes():
        pass

    
    def get_lum_nue():
        pass
    
    
    def get_lum_anue():
        pass

    
    def get_lum_nux():
        pass

    
    def get_lums():
        pass

    
    def get_aveene_nue():
        pass
    
    
    def get_aveene_anue():
        pass

    
    def get_aveene_nux():
        pass

    
    def get_aveenes():
        pass

    
    def get_number_spectra_nue(self):
        """ returns the array of nue_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_fluxes_nue

    
    def get_number_spectra_anue(self):
        """ returns the array of anue_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_fluxes_anue

    
    def get_number_spectra_nux(self):
        """ returns the array of nux_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_fluxes_nux
    
    
    def get_energy_spectra_nue(self):
        """ returns the array of nue_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_lums_nue

    
    def get_energy_spectra_anue(self):
        """ returns the array of anue_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_lums_anue


    def get_energy_spectra_nux(self):
        """ returns the array of nux_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__osc_lums_nux