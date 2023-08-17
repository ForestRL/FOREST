import sys
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
sys.path.append("..")
from nu_osc_model import nu_osc_model
from nu_reaction import nu_reaction
import numpy as np
from numpy import random
import tools
from event_generator import event_generator

class analytic_generator(event_generator):
    def __init__(self, M:float=1.4, R:float=10.0, gbeta:float=3.0, etot:float=1.0e53,
                  volume:float=32.5, distance:float=10.0):
        self.__M = M
        self.__R = R
        self.__volume = volume
        self.__distance = distance
        self.__gbeta = gbeta
        self.__etot = etot

        self.__t0  = 210.e0*(M/1.4e0)**1.2/(R/10.e0)**1.2*(gbeta/3.e0)**0.8/(etot/1.e52)**0.2
        self.__enes = np.linspace(0.01, 100, 200)


    def __get_eventrate_t(self, t:float):
   #     luminosity = 3.3e51*(self.__M/1.4)**(6)*(self.__R/10.0)**(-6)\
   #                       *(self.__gbeta/3.0)**(4.)*((t + self.__t0)/100.)**(-6.)
   #     ave_ene =  16*(self.__M/1.4)**(3./2)*(self.__R/10)**(-2.)*(self.__gbeta/3.)\
   #                     *((t+self.__t0)/100.)**(-3./2)

        rate = 720.*(self.__volume/32.5)*(self.__distance/10.)**(-2.)*(self.__R/10.)\
                *(self.__gbeta/3.)**(5.)*((t+self.__t0)/100.)**(-15./2)


    def get_eventrate(self, react_name:str):
        """
            returns eventrate of the reaction
                retrun array of eventrate
        """

        #luminosity = 3.3e51*(self.__M/1.4)**(6)*(self.__R/10.0)**(-6)\
        #                *(self.__gbeta/3.0)**(4)*(self.__t0 + )**(-5./2)
        ave_ene =  16*(self.__M/1.4)**(3./2)*(self.__R/10)**(-2.)*(self.__gbeta/3.)


    def get_tot_events(self, react_name:str, flavor:int, t_end:float):
        """
            returns total of the reaction
                retrun array of total events
        """
        pass


    def get_events(self, t_start:float, t_end:float, react_names=[]):
        """
            returns arrays of individual events
                retrun array of events
        """
        tot_start = 720.0*100*(-2.0/13)*(self.__volume/32.5)*(self.__distance/10.0)**(-2.0)*(self.__R/10.0)\
                *(self.__gbeta/3.0)**(5.0)*((t_start+self.__t0)/100.0)**(-13.0/2)

        tot_end = 720.0*100*(-2.0/13)*(self.__volume/32.5)*(self.__distance/10.0)**(-2.0)*(self.__R/10.0)\
                *(self.__gbeta/3.0)**(5.0)*((t_end+self.__t0)/100.0)**(-13.0/2)

        t_center = (t_start + t_end)/2.0
        kT = 5.0e0*(self.__M/1.4e0)**1.5/(self.__R/10.e0)**2*(self.__gbeta/3.e0)/((t_center+self.__t0)/100.e0)**(1.5e0)
        dn_exp = tot_end - tot_start
        dn = np.random.poisson(dn_exp, 1)
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []
        for time in times:
            F = tools.fermi(kT, 4, self.__enes)
            ene = tools.get_value_random(self.__enes, F)
            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)
            event_list.append({"time":time, "nu_ene":ene+1.3, "ev_ene":ene, "phi":phi, "theta":theta, "id":-11})
        return event_list


if __name__ == "__main__":

    gen = analytic_generator()

    times = np.linspace(0.01, 100, 10000)
    for time2, time1 in zip(times[1:], times[:-1]):
        gen.get_events(time1, time2)
