from abc import ABCMeta, abstractmethod, ABC
from spectra import spectra
import numpy as np

class analytic_spectra(spectra):
    """ addresses neutrino spectra
    """

    def __init__(self, M:float=1.4, R:float=10.0, gbeta:float=3.0, etot:float=1.0e53,
                 distance:float=8.0, tend:float=100.0, e_min:float=1.0, e_max:float=300, ebins:int=18) -> None:
        super().__init__()

        mev_to_erg = 1.60218e-6

        self.__distance = distance
        self.__times = np.arange(0.0, tend, 0.1)

        t0 = 210.e0*(M/1.4e0)**1.2/(R/10.e0)**1.2*(gbeta/3.e0)**0.8/(etot/1.e52)**0.2
        kT = 5.0e0*(M/1.4e0)**1.5/(R/10.e0)**2*(gbeta/3.e0)/((self.__times+t0)/100.e0)**(1.5e0)
        lum_anue_t =  3.3e51*(M/1.4)**(6)*(R/10.0)**(-6)\
                              *(gbeta/3.0)**(4.)*((self.__times + t0)/100.)**(-6.)
        ave_ene_anue_t = 16.0e0*(M/1.4e0)**1.5/(R/10.e0)**-2*(gbeta/3.e0)/((self.__times+t0)/100.e0)**(-1.5e0)
        flux_anue_t = lum_anue_t/(ave_ene_anue_t*mev_to_erg)

        r = (e_max/e_min)**(1.0/(ebins - 1))
        self.__ene2_bins = np.tile(np.array([ e_min*r**i for i in range(ebins)]), (len(self.__times),1))
        self.__ene1_bins = np.tile(np.append(np.array([0.0]), self.__ene2_bins[1,:-1]), (len(self.__times),1))
        self.__ene_cen = (self.__ene1_bins + self.__ene2_bins)*0.5
        self.__fluxes_nue = np.zeros_like(self.__ene1_bins)
        self.__fluxes_anue = self.__fermi(kT, self.__ene1_bins, self.__ene2_bins, 2)*np.tile(flux_anue_t, (ebins,1)).T
        self.__fluxes_nux = np.zeros_like(self.__ene1_bins)
        self.__lums_nue = np.zeros_like(self.__ene1_bins)
        self.__lums_anue = self.__fermi(kT, self.__ene1_bins, self.__ene2_bins, 3)*np.tile(lum_anue_t, (ebins,1)).T
        self.__lums_nux = np.zeros_like(self.__ene1_bins)


    def __fermi(self, kT, ene1_bins, ene2_bins, n):
        f = np.zeros((len(kT), len(ene1_bins[0,:])))
        ene_cen = (ene1_bins + ene2_bins)*0.5
        width = ene2_bins - ene1_bins
        kT2d = np.tile(kT, (len(ene1_bins[0,:]),1) ).T
        f = ene_cen**n/(np.exp(ene_cen/kT2d) + 1.0)
        F = np.tile(np.sum(f*width, axis=1), (len(ene1_bins[0,:]),1)).T
        return f/F
        


        


    def get_luminosities(self):
        """calculates luminosities for each flavor.
            return:
                numpy arrays
                time, nue lum, anue lum, nux_lum
        """
        pass

    
    def get_fluxes(self):
        """calculates fluxes for each flavor.
            return:
                numpy arrays
                time, nue flux, anue flux, nux_flux
        """
        pass 

    
    def get_average_energies(self):
        """ calculates average energies for each flavor.
            the formula in the Nakazato format is
            Luminsity/Flux/1.6022e-6
            return:
                numpy arrays
                time, nue ene, anue ene, nux ene
        """
        pass

    
    def get_ene1_bins(self):
        """ returns the array of ene1_bins.
            return:
                   numpy array
                   ene1_bins
        """
        return self.__ene1_bins

    
    def get_ene2_bins(self):
        """ returns the array of ene2_bins.
            return:
                   numpy array
                   ene2_bins
        """
        return self.__ene2_bins

    
    def get_ene_center_bins(self):
        """ returns the array of ene_center_bins.
            return:
                   numpy array
                   ene_center_bins
        """
        return self.__ene_cen
    
    
    def get_nue_energy_spectra(self):
        """ returns the array of nue_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__lums_nue


    
    def get_anue_energy_spectra(self):
        """ returns the array of anue_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__lums_anue

        
    def get_nux_energy_spectra(self):
        """ returns the array of nux_energy_spectra.
            return:
                   numpy array
                   
        """
        return self.__lums_nux

    
    def get_nue_number_spectra(self):
        """ returns the array of nue_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__fluxes_nue
    
    
    def get_anue_number_spectra(self):
        """ returns the array of anue_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__fluxes_anue

    
    def get_nux_number_spectra(self):
        """ returns the array of nux_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__fluxes_nux

    
    def get_times(self):
        """ returns the array of times.
            return:
                   nunpy array
        """
        return self.__times

    
    def get_distance(self) -> float:
        """ returns the distance of the SN from earth.
            return:
                   distance (kpc)
        """
        return self.__distance
