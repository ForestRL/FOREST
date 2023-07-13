from abc import ABCMeta, abstractmethod
import sys
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
sys.path.append("..")
from nu_osc_model import nu_osc_model
from nu_reaction import nu_reaction

class detector(metaclass=ABCMeta):
    """Abstract class for detectors"""

    @abstractmethod
    def __init__(self, nu_osc:nu_osc_model, nu_reactions:list[nu_reaction]):
        pass


    @abstractmethod
    def set_nu_reaction(self):
        pass


    @abstractmethod
    def set_nu_oscillation(self):
        pass


    @abstractmethod
    def set_background(self):
        pass


    @abstractmethod
    def get_event_rate(self) -> float:
        pass


    @abstractmethod
    def get_events(self) -> list[float]:
        pass

