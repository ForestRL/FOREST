from abc import ABCMeta, abstractmethod, ABC

class nu_reaction(metaclass=ABCMeta):
    """Abstract class for calculating neutrino reactions"""
    @abstractmethod
    def dcs_nue(self, dsigmadE, Enu, E_target) -> float:
        pass

    @abstractmethod
    def dcs_cos_nue(self, dcs, Enu, cos) -> float:
        pass

    @abstractmethod
    def cs_nue(self, dcs, Enu, cos) -> float:
        pass

    @abstractmethod
    def dcs_anue(self, dsigmadE, Enu, E_target) -> float:
        pass

    @abstractmethod
    def dcs_cos_anue(self, dcs, Enu, cos) -> float:
        pass

    @abstractmethod
    def cs_anue(self, dcs, Enu, cos) -> float:
        pass

    @abstractmethod
    def dcs_nux(self, dsigmadE, Enu, E_target) -> float:
        pass

    @abstractmethod
    def dcs_cos_nux(self, dcs, Enu, cos) -> float:
        pass

    @abstractmethod
    def cs_nux(self, dcs, Enu, cos) -> float:
        pass
    