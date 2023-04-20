from abc import ABCMeta, abstractmethod, ABC
import numpy as np

class nu_reaction(metaclass=ABCMeta):
    """Abstract class for calculating neutrino reactions"""
    @abstractmethod
    def dcs_nue(self, dsigmadE:float, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_nue(self, dcs:float, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_nue(self, dcs:float, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def dcs_anue(self, dsigmadE:float, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_anue(self, dcs:float, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_anue(self, dcs:float, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def dcs_nux(self, dsigmadE:float, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_nux(self, dcs:float, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_nux(self, dcs:float, Enu:float, cos:float) -> float:
        pass
    
    @abstractmethod
    def cs_table_anue(self, n_cos:int=200, n_e:int=300, E_min:float=0.0, E_max:float=300) -> list[np.array]:
        pass