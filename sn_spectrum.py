"""
sn_spctrum.py
Copyright (c) 2022 Masamitsu Mori

The sn_snpectrum module addresses a SN spectrum file written in the Nakazato format.

Nakazato format:
 "T0
  E0  E1    dN_nue(T0)/dE1  dN_anue(T0)/dE1  dN_nux(T0)/dE1  dL_nue(T0)/dE1  dL_anue(T0)/dE1  dL_nux(T0)/dE1 
  E1  E2    dN_nue(T0)/dE2  dN_anue(T0)/dE1  dN_nux(T0)/dE2  dL_nue(T0)/dE2  dL_anue(T0)/dE2  dL_nux(T0)/dE2 
  .
  .
  .
  E19 E20   dN_nue(T0)/dE20 dN_anue(T0)/dE20 dN_nux(T0)/dE20 dL_nue(T0)/dE20 dL_anue(T0)/dE20 dL_nux(T0)/dE20

  T1
  E0  E1    dN_nue(T0)/dE1  dN_anue(T1)/dE1  dN_nux(T1)/dE1  dL_nue(T1)/dE1  dL_anue(T1)/dE1  dL_nux(T1)/dE1
  .
  .
  .
",where Tn [s] is a time measured from the bounce and Ek [MeV] is a neutrino energy.

 When launched in the command line, the module calculates luminosity and average energy.
"""

import numpy as np
import sys


class SNspectrum:
    """ addresses neutrino spectra
    """


    def get_luminosities(self):
        """calculates luminosities for each flavor.
            return:
                numpy arrays
                time, nue lum, anue lum, nux_lum
        """

        # calculate flux(t) = sum_k=1,Nbins((E_k - E_{k-1})*(dN/dE)) for each time
        lum_nue_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__lums_nue, axis=1)
        lum_anue_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__lums_anue, axis=1)
        lum_nux_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__lums_nux, axis=1)

        # These sentences correspond to below:
        #lum_nue_array = np.empty((0))
        #lum_anue_array = np.empty((0))
        #lum_nux_array = np.empty((0))
        #for e1, e2, dL_dedt_nue, dL_dedt_anue, dL_dedt_nux in zip(self.__ene1_bins[:], self.__ene2_bins[:], self.__lums_nue[:], self.__lums_anue[:], self.__lums_nux[:]):
            
        #    lum_nue = np.sum((e2-e1)*dL_dedt_nue) 
        #    lum_anue = np.sum((e2-e1)*dL_dedt_anue) 
        #    lum_nux = np.sum((e2-e1)*dL_dedt_nux) 

        #    lum_nue_array = np.hstack((lum_nue_array, lum_nue))
        #    lum_anue_array = np.hstack((lum_anue_array, lum_anue))
        #    lum_nux_array = np.hstack((lum_nux_array, lum_nux))

        return self.__times, lum_nue_array, lum_anue_array, lum_nux_array
 

    def get_fluxes(self):
        """calculates fluxes for each flavor.
            return:
                numpy arrays
                time, nue flux, anue flux, nux_flux
        """

        # calculate flux(t) = sum_k=1,Nbins((E_k - E_{k-1})*(dN/dE)) for each time
        flux_nue_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__fluxes_nue, axis=1)
        flux_anue_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__fluxes_anue, axis=1)
        flux_nux_array = np.sum((self.__ene2_bins-self.__ene1_bins)*self.__fluxes_nux, axis=1)

        # These sentences correspond to below:
        #flux_nue_array = np.empty((0))
        #flux_anue_array = np.empty((0))
        #flux_nux_array = np.empty((0))
        #for e1, e2, dN_dedt_nue, dN_dedt_anue, dN_dedt_nux in zip(self.__ene1_bins[:], self.__ene2_bins[:], self.__fluxes_nue[:], self.__fluxes_anue[:], self.__fluxes_nux[:]):
            
        #    flux_nue = np.sum((e2-e1)*dN_dedt_nue) 
        #    flux_anue = np.sum((e2-e1)*dN_dedt_anue) 
        #    flux_nux = np.sum((e2-e1)*dN_dedt_nux) 

        #   flux_nue_array = np.hstack((flux_nue_array, flux_nue))
        #   flux_anue_array = np.hstack((flux_anue_array, flux_anue))
        #   flux_nux_array = np.hstack((flux_nux_array, flux_nux))

        return self.__times, flux_nue_array, flux_anue_array, flux_nux_array
    

    def get_average_energies(self):
        """ calculates average energies for each flavor.
            the formula in the Nakazato format is
            Luminsity/Flux/1.6022e-6
            return:
                numpy arrays
                time, nue ene, anue ene, nux ene
        """
        erg_to_MeV = 1./1.6022e-6
        times, flux_nue, flux_anue, flux_nux = self.get_fluxes()
        times, lum_nue, lum_anue, lum_nux    = self.get_luminosities()

        return times, lum_nue/flux_nue*erg_to_MeV,  lum_anue/flux_anue*erg_to_MeV,  lum_nux/flux_nux*erg_to_MeV

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
        return (self.__ene2_bins + self.__ene1_bins)/2.0


    def get_anue_number_spectra(self):
        """ returns the array of anue_number_spectra.
            return:
                   numpy array
                   
        """
        return self.__fluxes_anue


    def get_times(self):
        """ returns the array of times.
            return:
                   nunpy array
        """

        return self.__times


    def  __init__(self, spectrum_file: str):
        """initializes the sn_spectrum class.
            loads sn spectrum file.

            arg:
                spectrum_file path
        """

        self.__times = np.empty((0)) # array of Tn
        self.__n_ene_bin = 20 # number of energy bins
        self.__ene1_bins = np.empty((0, self.__n_ene_bin)) # array of first energy bins
        self.__ene2_bins = np.empty((0, self.__n_ene_bin)) # array of second energy bins
        self.__fluxes_nue   = np.empty((0, self.__n_ene_bin)) # array of flux of nue
        self.__fluxes_anue  = np.empty((0, self.__n_ene_bin)) # array of flux of anue
        self.__fluxes_nux   = np.empty((0, self.__n_ene_bin)) # array of flux of nux
        self.__lums_nue     = np.empty((0, self.__n_ene_bin)) # array of luminosity of nue
        self.__lums_anue    = np.empty((0, self.__n_ene_bin)) # array of luminosity of anue
        self.__lums_nux     = np.empty((0, self.__n_ene_bin)) # array of luminosity of nux

        with open(spectrum_file, mode='r') as f:
            
            while True:
                line = f.readline()
                if line == '':
                    break

                time = float(line)
                self.__times = np.hstack((self.__times, np.array((time))))
                
                ene1_bin = []
                ene2_bin = []
                flux_nue  = []
                flux_anue = []
                flux_nux  = []
                lum_nue   = []
                lum_anue  = []
                lum_nux   = []

                for i in range(self.__n_ene_bin):
                    data = f.readline().split()
                    
                    ene1_bin.append(float(data[0]))
                    ene2_bin.append(float(data[1]))
                    flux_nue.append(float(data[2]))
                    flux_anue.append(float(data[3]))
                    flux_nux.append(float(data[4]))
                    lum_nue.append(float(data[5]))
                    lum_anue.append(float(data[6]))
                    lum_nux.append(float(data[7]))

                self.__ene1_bins = np.vstack((self.__ene1_bins, np.array(ene1_bin)))
                self.__ene2_bins = np.vstack((self.__ene2_bins, np.array(ene2_bin)))
                self.__fluxes_nue = np.vstack((self.__fluxes_nue, np.array(flux_nue)))
                self.__fluxes_anue = np.vstack((self.__fluxes_anue, np.array(flux_anue)))
                self.__fluxes_nux = np.vstack((self.__fluxes_nux, np.array(flux_nux)))
                self.__lums_nue = np.vstack((self.__lums_nue, np.array(lum_nue)))
                self.__lums_anue = np.vstack((self.__lums_anue, np.array(lum_anue)))
                self.__lums_nux = np.vstack((self.__lums_nux, np.array(lum_nux)))

                f.readline()


# If launched in command line
if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print('Usage: python3 ' + sys.argv[0] + ' ' + '[SN spectrum file]')
        print('Output: lum.dat, ave.dat, flux.dat')
        print('Output format:\n #1 time #2 electron neutrino   #3 anti-electron neutrino  #4 other neutrinos')
        sys.exit(-1)

    sn_spectrum = SNspectrum(sys.argv[1])

    with open('lum.dat', 'w') as f:
        time_array, lum_nue_array, lum_anue_array, lum_nux_array = sn_spectrum.get_luminosities()
        for time, lum_nue, lum_anue, lum_nux in zip(time_array, lum_nue_array, lum_anue_array, lum_nux_array):
            f.write(str(time) + ' ' + str(lum_nue) + ' ' + str(lum_anue) + ' ' + str(lum_nux))

    with open('flux.dat', 'w') as f:
        time_array, flux_nue_array, flux_anue_array, flux_nux_array = sn_spectrum.get_average_energies()
        for time, flux_nue, flux_anue, flux_nux in zip(time_array, flux_nue_array, flux_anue_array, flux_nux_array):
            f.write(str(time) + ' ' + str(flux_nue) + ' ' + str(flux_anue) + ' ' + str(flux_nux))
    
    with open('ave.dat', 'w') as f:
        time_array, ave_nue_array, ave_anue_array, ave_nux_array = sn_spectrum.get_average_energies()
        for time, ave_nue, ave_anue, ave_nux in zip(time_array, ave_nue_array, ave_anue_array, ave_nux_array):
            f.write(str(time) + ' ' + str(ave_nue) + ' ' + str(ave_anue) + ' ' + str(ave_nux))
