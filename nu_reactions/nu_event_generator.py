import sys
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
sys.path.append("..")
from nu_osc_model import nu_osc_model
from nu_reaction import nu_reaction
import numpy as np
from numpy import random
import tools

class nu_event_generator:
    def __init__(self, nu_osc:nu_osc_model, nu_reaction_list:list[nu_reaction], nu_target_list:list[float], n_cos_res:int):

        self.__kpc_cm = 3.086e21 # 1 kpc in cm

        self.__nu_osc_model = nu_osc
        self.__sn_spectra = nu_osc.get_sn_spectra()
        self.__nu_reaction_dict = {}
        self.__nu_target_dict = {}


        for nu_react, nu_target in zip(nu_reaction_list, nu_target_list):
            self.__nu_reaction_dict[nu_react.get_reaction_name()] = nu_react
            self.__nu_target_dict[nu_react.get_reaction_name()] = nu_target
        

        self.__ene_centers_all = nu_osc.get_sn_spectra().get_ene_center_bins()
        self.__cs_ene = {}
        self.__dcsdEdcos = {}

        self.__ene1_bins_all = self.__sn_spectra.get_ene1_bins()
        self.__ene2_bins_all = self.__sn_spectra.get_ene2_bins()
        self.__ene_bins_all = np.vstack((self.__ene1_bins_all[:,0], self.__ene2_bins_all.T)).T
        self.__ene_widths = self.__ene2_bins_all - self.__ene1_bins_all

        cos_list = np.linspace(-1.0, 1.0, n_cos_res)
        self.__cos_list = cos_list

        for nu_react in nu_reaction_list:
            self.__cs_ene[nu_react.get_reaction_name()] = [nu_react.cs_nue(self.__ene_bins_all),
                                                            nu_react.cs_anue(self.__ene_bins_all),
                                                            nu_react.cs_nux(self.__ene_bins_all)]
            self.__dcsdEdcos[nu_react.get_reaction_name()] = [nu_react.dcs_cos_nue(self.__ene_bins_all, cos_list).transpose(1,2,0),
                                                              nu_react.dcs_cos_anue(self.__ene_bins_all, cos_list).transpose(1,2,0),
                                                              nu_react.dcs_cos_nux(self.__ene_bins_all, cos_list).transpose(1,2,0)]


        self.__event_spectra = {}
        #event_spectra_temp = np.empty((len(self.__ene_centers_all[:,0]), len(self.__ene_bins_all[0,:]), 3), dtype=float)
        for nu_react in nu_reaction_list:

            spectra_nue_temp = nu_osc.get_number_spectra_nue()*self.__cs_ene[nu_react.get_reaction_name()][0][:,1:]
            spectra_anue_temp = nu_osc.get_number_spectra_anue()*self.__cs_ene[nu_react.get_reaction_name()][1][:,1:]
            spectra_nux_temp = nu_osc.get_number_spectra_nux()*self.__cs_ene[nu_react.get_reaction_name()][2][:,1:]

            event_spectra_temp = [np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_nue_temp.T)).T
                                         ,np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_anue_temp.T)).T
                                         ,np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_nux_temp.T)).T]


            #event_spectra_temp[:,:,0] = np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_nue_temp.T)).T
            #event_spectra_temp[:,:,1] = np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_anue_temp.T)).T
            #event_spectra_temp[:,:,2] = np.vstack((np.zeros(len(self.__ene_bins_all[:,0])), spectra_nux_temp.T)).T
            
            self.__event_spectra[nu_react.get_reaction_name()] = event_spectra_temp

        self.__event_rates = {}
        self.__cum_events = {}

        for nu_react in nu_reaction_list:
            
            dN_dE_nue = nu_osc.get_number_spectra_anue()*self.__nu_target_dict[nu_react.get_reaction_name()]/(4.0*np.pi*self.__sn_spectra.get_distance()**2*self.__kpc_cm**2)\
                                            *(self.__cs_ene[nu_react.get_reaction_name()][0][:,1:]+ self.__cs_ene[nu_react.get_reaction_name()][0][:,:-1])*0.5*self.__ene_widths
            dN_dE_anue = nu_osc.get_number_spectra_anue()*self.__nu_target_dict[nu_react.get_reaction_name()]/(4.0*np.pi*self.__sn_spectra.get_distance()**2*self.__kpc_cm**2)\
                                            *(self.__cs_ene[nu_react.get_reaction_name()][1][:,1:]+ self.__cs_ene[nu_react.get_reaction_name()][1][:,:-1])*0.5*self.__ene_widths
            dN_dE_nux = nu_osc.get_number_spectra_anue()*self.__nu_target_dict[nu_react.get_reaction_name()]/(4.0*np.pi*self.__sn_spectra.get_distance()**2*self.__kpc_cm**2)\
                                            *(self.__cs_ene[nu_react.get_reaction_name()][2][:,1:]+ self.__cs_ene[nu_react.get_reaction_name()][2][:,:-1])*0.5*self.__ene_widths
            
            self.__event_rates[nu_react.get_reaction_name()] = [np.sum(dN_dE_nue, axis=1),
                                                                np.sum(dN_dE_anue, axis=1),
                                                                np.sum(dN_dE_nux, axis=1)]
#            self.__cum_events[nu_react.get_reaction_name()] = self.__event_rates[nu_react.get_reaction_name()].cumsum()


    def get_eventrate(self, react_name:str):
        """
            returns eventrate of the reaction
                retrun array of eventrate
        """

        return self.__event_rates[react_name]


    def get_tot_events(self, react_name:str, flavor:int, t_end:float):
        time = self.__sn_spectra.get_times()

        rate = self.__event_rates[react_name][flavor]
        cum_events = np.append(0.0, np.cumsum((rate[1:] + rate[:-1])*(time[1:] - time[:-1])*0.5))

        if(t_end >= time[-1]):
            return cum_events[-1]

        max_i = np.argmax(time[time - t_end <= 0])

        time_1 = time[max_i]
        time_2 = time[max_i + 1]

        return (cum_events[max_i+1] - cum_events[max_i])*(time_2 - time_1)*(t_end - time_1) \
                + cum_events[max_i]


    def get_events(self, t_start:float, t_end:float, react_names:list[str]=[]):
        
        for i in range(3):
            tot_start = self.get_tot_events(react_names[0], i, t_start)
            tot_end = self.get_tot_events(react_names[0], i, t_end)
            dn_exp = tot_end - tot_start
 
            dn = random.poisson(dn_exp, 1)
            
            times = np.random.rand(dn[0])*(t_end - t_start) + t_start
            times.sort()
            for t in times:

                event_spectrum_t = tools.interpolate_2d_array(t, self.__sn_spectra.get_times(), self.__event_spectra[react_names[0]][i])
                ene_bins_t = tools.interpolate_2d_array(t, self.__sn_spectra.get_times(), self.__ene_bins_all)
                dcsdEdcos_t = tools.interpolate_3d_array(t, self.__sn_spectra.get_times(), self.__dcsdEdcos[react_names[0]][i])

                ene = tools.get_value_random(ene_bins_t, event_spectrum_t)
                dcs_dcos = tools.interpolate_2d_array(ene, ene_bins_t, dcsdEdcos_t)
                cos = tools.get_value_random(self.__cos_list , dcs_dcos)

                pid = self.__nu_reaction_dict[react_names[0]].get_PID1()
                Ee = self.__nu_reaction_dict[react_names[0]].get_outgoing_particle_energy(cos, ene)

                print(t, ene, cos, Ee, pid)


    def __init_cs_table(self):
        pass


if __name__ == "__main__":
    sys.path.append("../")
    sys.path.append("../nu_osc_models")
    from sn_spectra import SNspectra
    from MSW_normal import MSW_normal
    from inverse_beta_decay import inverse_beta_decay

    spectra = SNspectra("../z9.6_ver2_nuspectra_nonzero.dat")
    nu_osc = MSW_normal(spectra)
    ibd = inverse_beta_decay()

    ev_gen = nu_event_generator(nu_osc, [ibd], [2.173e33], 200)
    times = spectra.get_times()
    for time1, rate in zip(times[1:], ev_gen.get_eventrate(ibd.get_reaction_name())[1][1:]):
        print(time1, rate)

