from detector import detector
import numpy as np
import sys
sys.path.append("../")
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
from nu_osc_model import nu_osc_model
from nu_reactions.event_generator import event_generator
from nu_reactions.analytic_generator import analytic_generator
import tools
import typing

class super_kamiokande(detector):
    """Impelement for Super-Kamiokande"""

    def __init__(self, nu_osc:nu_osc_model, ev_gen:event_generator):
        self.__height = 41.4
        self.__diameter = 39.3
        self.__fv_distance = 2.0
        self.__n_protons = 2.173e33
        self.__env_slope = 10
        self.__env_center = 39.3

        self.__full_volume = np.pi*(self.__diameter/2.0)**2*self.__height
        self.__top_bottom_volume = np.pi*(self.__diameter/2.0)**2*(2*self.__fv_distance)
        self.__cylinder_volume = self.__full_volume - np.pi*(self.__diameter/2.0 - self.__fv_distance)**2*self.__height
        self.__fiducial_volume = np.pi*(self.__diameter/2.0 - self.__fv_distance)**2*(self.__height - 2*self.__fv_distance)

        self.__ev_gen = ev_gen

        self.__bg_data = np.loadtxt("sk_bg.dat")

        self.__bg_data[:,1:] = np.where(self.__bg_data[:,1:] > 1e-100, self.__bg_data[:,1:], 1e-100)
#        self.__bg_data[:,3] = np.where(self.__bg_data[:,3] > self.__bg_data[:,1], self.__bg_data[:,1], self.__bg_data[:,3])
        self.__bin_width = self.__bg_data[1,0] - self.__bg_data[0,0] 
        self.__bin_lowedges = self.__bg_data[:,0] - self.__bin_width
        self.__bin_highedges = self.__bg_data[:,0] + self.__bin_width

        self.__r2_outFV = np.linspace((self.__diameter - self.__fv_distance)**2, self.__diameter**2, 200)
        self.__z_outFV_p = np.linspace( self.__height/2.0-self.__fv_distance, self.__height/2.0, 200)

        self.__make_envelope()


    def __make_envelope(self):

        r_fv2 = (self.__diameter - self.__fv_distance)**2
        r_wall2 = self.__diameter**2
        s2 = self.__env_slope**2
        top_wall = self.__height/2.0
        top_fv = self.__height/2.0 - self.__fv_distance

        ar2 = 1.0/(s2*(np.exp(r_wall2/s2)-np.exp(r_fv2/s2)))
        self.__ar2 = ar2

        az = 1.0/(self.__env_slope*(np.exp(top_wall/self.__env_slope)-np.exp(top_fv/self.__env_slope)))
        self.__az = az

    def __bg_model_func(self, a, x):
        return a*np.exp(x/self.__env_slope**2)


    def generate_bg_in_fv(self, t_start:float, t_end:float):
        tot_rate = np.sum(self.__bg_data[:,1]*self.__bin_width)*self.__fiducial_volume
        dn_exp = (t_end - t_start)*tot_rate
        dn = np.random.poisson(dn_exp, 1)
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []

        for t in times:
            theta = 2*np.pi*np.random.rand()
            cos = np.cos(theta)
            sin = np.sin(theta)
            r = np.sqrt((self.__diameter/2 - self.__fv_distance)**2*np.random.rand())
            x = r*cos
            y = r*sin
            z = (np.random.rand() - 0.5)*self.__height
            
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,1])
            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)        
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0})
            #print(t, ene, r*cos, r*sin, z)
        return event_list


    def generate_bg_out_fv_top_bottom(self, t_start:float, t_end:float):
        tot_rate = np.sum(self.__bg_data[:,3]*self.__bin_width)*self.__top_bottom_volume
        dn_exp = (t_end - t_start)*tot_rate
        dn = np.random.poisson(dn_exp, 1)
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []

        for t in times:
            theta = 2*np.pi*np.random.rand()
            cos = np.cos(theta)
            sin = np.sin(theta)
            r2 = np.random.rand()*self.__diameter**2
            
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,1])
            
            pdf = self.__az*np.exp(self.__z_outFV_p/self.__env_slope)
            z = tools.get_value_random(self.__r2_outFV, pdf)
            z *= np.random.choice([1,-1])

            x = np.sqrt(r2)*cos
            y = np.sqrt(r2)*sin

            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)        
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0})
            #print(t, ene, x, y, z)

        return event_list


    def generate_bg_out_fv_cylinder(self, t_start:float, t_end:float):
        tot_rate = np.sum(self.__bg_data[:,3]*self.__bin_width)*self.__cylinder_volume
        dn_exp = (t_end - t_start)*tot_rate
        dn = np.random.poisson(dn_exp, 1)
        
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []
        for t in times:
            theta = 2*np.pi*np.random.rand()
            cos = np.cos(theta)
            sin = np.sin(theta)
            z = (np.random.rand() - 0.5)*self.__height
            
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,1])
            
            #a = tools.interpolate_1d_array(ene, self.__bg_data[:,0], self.__a)
            pdf = self.__ar2*np.exp(self.__r2_outFV/self.__env_slope**2)
            r2 = tools.get_value_random(self.__r2_outFV, pdf)

            x = np.sqrt(r2)*cos
            y = np.sqrt(r2)*sin
            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0})
            #print(t, ene, x, y, z)

        return event_list


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


    def get_events(self, t_start:float, t_end:float) -> typing.List[float]:
        bg_in_fv = self.generate_bg_in_fv(t_start, t_end)
        bg_out_fv_cylinder = self.generate_bg_out_fv_cylinder(t_start, t_end)
        bg_out_fv_top_bottom = self.generate_bg_out_fv_top_bottom(t_start, t_end)
        sn_ev = self.__ev_gen.get_events(t_start, t_end)

        ev_list = []
        ev_list.extend(bg_in_fv)
        ev_list.extend(bg_out_fv_cylinder)
        ev_list.extend(bg_out_fv_top_bottom)
        ev_list.extend(sn_ev)
        for ev in ev_list:
            print(ev["time"], ev["ev_ene"], ev["id"])

if __name__ == "__main__":
    import sys
    sys.path.append("../")
    sys.path.append("../nu_osc_models")
    sys.path.append("../nu_reactions")
    from sn_spectra import SNspectra
    from MSW_normal import MSW_normal
    from inverse_beta_decay import inverse_beta_decay
    from nu_event_generator import nu_event_generator


 #   ev_gen = nu_event_generator(nu_osc, [ibd], [2.173e33], 200)
    ev_gen_ana = analytic_generator(M=2.0)
    sk = super_kamiokande(None, ev_gen_ana)
    times = np.linspace(0.01, 200, 20000)
    for time1, time2 in zip(times[1:], times[:-1]):
        sk.get_events(time2, time1)
