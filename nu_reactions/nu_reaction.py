from abc import ABCMeta, abstractmethod, ABC
import numpy as np

class nu_reaction(metaclass=ABCMeta):
    """Abstract class for calculating neutrino reactions"""

    @abstractmethod
    def get_reaction_name(self) -> str:
        pass

    @abstractmethod
    def dcs_nue(self, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_nue(self, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_nue(self, Enu:float) -> float:
        pass

    @abstractmethod
    def dcs_anue(self, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_anue(self, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_anue(self, Enu:float) -> float:
        pass

    @abstractmethod
    def dcs_nux(self, Enu:float, E_target:float) -> float:
        pass

    @abstractmethod
    def dcs_cos_nux(self, Enu:float, cos:float) -> float:
        pass

    @abstractmethod
    def cs_nux(self, Enu:float) -> float:
        pass
    
    @abstractmethod
    def cs_table_anue(self, n_cos:int=200, n_e:int=300, E_min:float=0.0, E_max:float=300) -> list[np.array]:
        pass