from abc import ABCMeta, abstractmethod
import sys
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
sys.path.append("..")
from nu_osc_models.nu_osc_model import nu_osc_model
from nu_reactions.event_generator import event_generator
import typing

class detector(metaclass=ABCMeta):
    """Abstract class for detectors"""

    @abstractmethod
    def __init__(self, nu_osc:nu_osc_model, ev_gen:event_generator, background:bool=True, ene_cut:float=0.0):
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
    def get_events(self) -> typing.List[float]:
        pass

