from abc import ABCMeta, abstractmethod, ABC
import numpy as np

class nu_reaction(metaclass=ABCMeta):
    """Abstract class for calculating neutrino reactions"""

    @abstractmethod
    def get_reaction_name(self) -> str:
        pass

    @abstractmethod
    def dcs_nue(self, Enu, E_target):
        pass

    @abstractmethod
    def dcs_cos_nue(self, Enu, cos):
        pass

    @abstractmethod
    def cs_nue(self, Enu):
        pass

    @abstractmethod
    def dcs_anue(self, Enu, E_target):
        pass

    @abstractmethod
    def dcs_cos_anue(self, Enu, cos):
        pass

    @abstractmethod
    def cs_anue(self, Enu):
        pass

    @abstractmethod
    def dcs_nux(self, Enu, E_target):
        pass

    @abstractmethod
    def dcs_cos_nux(self, Enu, cos):
        pass

    @abstractmethod
    def cs_nux(self, Enu):
        pass
    
    @abstractmethod
    def cs_table_anue(self, n_cos:int=200, n_e:int=300, E_min:float=0.0, E_max:float=300) -> list[np.array]:
        pass

    @abstractmethod
    def get_PID1(self) ->int:
        pass

    @abstractmethod
    def get_outgoing_particle_energy(self, cos:float, Enu:float) -> float:
        pass