from detector import detector
import numpy as np
import sys
sys.path.append("../")
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
from nu_osc_model import nu_osc_model
from nu_reaction import nu_reaction

class super_kamiokande(detector):
    """Impelement for Super-Kamiokande"""

    def __init__(self, nu_osc:nu_osc_model, nu_reactions:list[nu_reaction]):
        pass

    def set_nu_reaction(self):
        pass


    def set_nu_oscillation(self):
        pass

    def set_background(self):
        bg_data = np.loadtxt("sk_bg.dat")
        bin_width = bg_data[1,0] - bg_data[0,0] 
        
        

        pass


    def get_event_rate(self) -> float:
        pass


    def get_events(self) -> list[float]:
        pass


if __name__ == "__main__":
    import sys
    sys.path.append("../")
    sys.path.append("../nu_osc_models")
    sys.path.append("../nu_reactions")
    from sn_spectra import SNspectra
    from MSW_normal import MSW_normal
    from inverse_beta_decay import inverse_beta_decay
    from nu_event_generator import nu_event_generator

    spectra = SNspectra("../z9.6_ver2_nuspectra_nonzero.dat")
    nu_osc = MSW_normal(spectra)
    ibd = inverse_beta_decay()

    ev_gen = nu_event_generator(nu_osc, [ibd], [2.173e33], 200)
    times = spectra.get_times()

    for time1, time2 in zip(times[1:], times[:-1]):
        ev_gen.get_events(time2, time1, [ibd.get_reaction_name()])
