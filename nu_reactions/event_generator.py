from abc import ABCMeta, abstractmethod, ABC

class event_generator(metaclass=ABCMeta):

    @abstractmethod
    def get_eventrate(self, react_name:str):
        """
            returns eventrate of the reaction
                retrun array of eventrate
        """
        pass

    @abstractmethod
    def get_tot_events(self, react_name:str, flavor:int, t_end:float):
        """
            returns total of the reaction
                retrun array of total events
        """
        pass

    @abstractmethod
    def get_events(self, t_start:float, t_end:float, react_names=[]):
        """
            returns arrays of individual events
                retrun array of events
        """

