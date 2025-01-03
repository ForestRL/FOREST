from detectors.detector import detector
import numpy as np
import sys
sys.path.append("../")
sys.path.append("../nu_osc_models")
sys.path.append("../nu_reactions")
from nu_osc_models.nu_osc_model import nu_osc_model
from nu_reactions.event_generator import event_generator
import tools
import typing

class super_kamiokande(detector):
    """Impelement for Super-Kamiokande"""
    VOLUME = 32.48 # SK inner volume in kton
    PROTONS =  2.173e33
    
    def __init__(self, nu_osc:nu_osc_model, ev_gen:event_generator, background:bool=True, ene_cut:float=4.5, ene_res:bool=False):
        self.__height = 36.2
        self.__diameter = 33.8
        self.__fv_distance = 2.0
        self.__n_protons = 2.173e33
        self.__env_slope = 5
        self.__env_center = 39.3
        
        if(ene_res):
            self.__ene_resolution = 1.5 # Mori et al. (2022)
        else:
            self.__ene_resolution = 0.0

        self.__full_volume = np.pi*(self.__diameter/2.0)**2*self.__height
        self.__top_bottom_volume = np.pi*(self.__diameter/2.0)**2*(2*self.__fv_distance)
        self.__cylinder_volume = self.__full_volume - np.pi*(self.__diameter/2.0 - self.__fv_distance)**2*self.__height
        self.__fiducial_volume = np.pi*(self.__diameter/2.0 - self.__fv_distance)**2*(self.__height - 2*self.__fv_distance)
#        self.__top_bottom_midle_volume = np.pi*(self.__diameter/2.0 - self.__fv_distance)**2*(2*self.__fv_distance)
#        self.__cylinder_midle_volume = np.pi*()

        self.__ev_gen = ev_gen

        self.__bg_data = np.loadtxt("sk_bg_extra_ver3.dat")

        self.__bg_data[:,1:] = np.where(self.__bg_data[:,1:] > 1e-100, self.__bg_data[:,1:], 1e-100)
        if(not background):
            self.__bg_data[:,1:] = 1e-100

        self.__bg_data = self.__bg_data[self.__bg_data[:,0] >= ene_cut]

        self.__sum_in = self.__bg_data[:,1].sum()
        self.__sum_out = self.__bg_data[:,3].sum()

        self.__bin_width = self.__bg_data[1,0] - self.__bg_data[0,0] 
        self.__bin_lowedges = self.__bg_data[:,0] - self.__bin_width
        self.__bin_highedges = self.__bg_data[:,0] + self.__bin_width

        self.__r2_outFV = np.linspace((self.__diameter/2.0 - self.__fv_distance)**2, (self.__diameter/2.0)**2, 200)
        self.__z_outFV_p = np.linspace( self.__height/2.0-self.__fv_distance, self.__height/2.0, 200)

        self.__make_envelope()
    

    def __make_envelope(self):
        r_fv2 = (self.__diameter/2.0 - self.__fv_distance)**2
        r_wall2 = self.__diameter**2
        a_cy = 1.0/5**2
        top_wall = self.__height/2.0
        top_fv = self.__height/2.0 - self.__fv_distance

        a_top = 1.5

        c_cy = a_cy/(np.exp(a_cy*r_wall2)-np.exp(a_cy*r_fv2))
        self.__pdf_slinder = c_cy*np.exp(self.__r2_outFV*a_cy)

        c_top = a_top/(np.exp(a_top*top_wall)-np.exp(a_top*top_fv))
        self.__pdf_top_bottom = c_top*np.exp(a_top*self.__z_outFV_p)


    def __make_envelope2(self):         
        c1 = self.sum_in/self.sum_out
        s2 = self.__env_slope**2
        r_fv2 = (self.__diameter/2.0 - self.__fv_distance)**2
        r_wall2 = self.__diameter**2
        top_wall = self.__height/2.0
        top_fv = self.__height/2.0 - self.__fv_distance

        f_int = s2*np.log(2*np.sinh((r_wall2-r_fv2)/s2)) 
        c0 = (1.0 - c1*(r_wall2 - r_fv2))/f_int 

        self.__pdf_slinder = c0*np.tanh((self.__r2_outFV - r_fv2)/s2)

        f_int = self.__env_slope*np.log(2*np.sinh((top_wall-top_fv)/self.__env_slope)) 
        c0 = (1.0 - c1*(top_wall - top_fv))/f_int 

        self.__pdf_top_bottom = c0*np.tanh((self.__z_outFV_p - top_fv)/self.__env_slope)
        
        
    def __bg_model_func(self, a, x):
        return a*np.exp(x/self.__env_slope**2)


    def __check_FV(self, x:float, y:float, z:float) -> bool:
        r = np.sqrt(x*x + y*y)

        if(r < self.__diameter/2 - self.__fv_distance and 
           np.abs(z) < self.__height/2 - self.__fv_distance):
            return True
        else:
            return False


    def generate_bg_in_fv(self, t_start:float, t_end:float):
        tot_rate = np.sum(self.__bg_data[:,1]*self.__bin_width)*self.__fiducial_volume
        dn_exp = (t_end - t_start)*tot_rate
        dn = np.random.poisson(dn_exp, 1)
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []

        for t in times:
            theta = 2*np.pi*np.random.rand() - np.pi
            cos = np.cos(theta)
            sin = np.sin(theta)
            r = np.sqrt((self.__diameter/2 - self.__fv_distance)**2*np.random.rand())
            x = r*cos
            y = r*sin
            z = (np.random.rand() - 0.5)*(self.__height - 2*self.__fv_distance)
            
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,1])
            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)        
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0, "fv":1})
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
            theta = 2*np.pi*np.random.rand() - np.pi
            cos = np.cos(theta)
            sin = np.sin(theta)
            r2 = np.random.rand()*(self.__diameter/2)**2
            
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,3])
            
            #pdf = self.__az*np.exp(self.__z_outFV_p/self.__env_slope)
            z = tools.get_value_random(self.__z_outFV_p, self.__pdf_top_bottom)
            z *= np.random.choice([1,-1])

            x = np.sqrt(r2)*cos
            y = np.sqrt(r2)*sin

            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)        
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0, "fv":0})
#            print(t, ene, x, y, z)

        return event_list


    def generate_bg_out_fv_cylinder(self, t_start:float, t_end:float):
        tot_rate = np.sum(self.__bg_data[:,3]*self.__bin_width)*self.__cylinder_volume
        dn_exp = (t_end - t_start)*tot_rate
        dn = np.random.poisson(dn_exp, 1)
        
        times = np.random.rand(dn[0])*(t_end - t_start) + t_start
        times.sort()
        event_list = []
        for t in times:
            theta = 2*np.pi*np.random.rand() - np.pi
            cos = np.cos(theta)
            sin = np.sin(theta)
            z = (np.random.rand() - 0.5)*self.__height
        
            ene = tools.get_value_random_hist(self.__bin_lowedges, self.__bin_highedges, self.__bg_data[:,3])

            #a = tools.interpolate_1d_array(ene, self.__bg_data[:,0], self.__a)
            #pdf = self.__ar2*np.exp(self.__r2_outFV/self.__env_slope**2)
            r2 = tools.get_value_random(self.__r2_outFV, self.__pdf_slinder)

            x = np.sqrt(r2)*cos
            y = np.sqrt(r2)*sin
            phi = 2.0*np.pi*np.random.rand()
            theta = np.arccos(-2.0*np.random.rand() + 1.0)
            event_list.append({"time":t, "ev_ene":ene, "nu_ene":0.0, "theta":theta, "phi":phi, "x":x, "y":y, "z":z, "id":0, "fv":0})
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
        sn_ev = self.__ev_gen.get_events(t_start, t_end, self.__ev_gen.get_react_names())
        
        for ev in sn_ev:
            theta = 2.0*np.pi*np.random.rand()
            r = np.sqrt(np.random.rand()*((self.__diameter/2.0)**2))
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            z = (np.random.rand() - 0.5)*self.__height
            ev["x"] = x
            ev["y"] = y
            ev["z"] = z

            ev["ev_ene"] = np.random.normal(ev["ev_ene"], self.__ene_resolution)
            
            if(ev["ev_ene"] < 0):
                ev["ev_ene"] = 0.0

            if(self.__check_FV(x,y,z)):
                ev["fv"] = 1
            else:
                ev["fv"] = 0


        ev_list = []
        ev_list.extend(bg_in_fv)
        ev_list.extend(bg_out_fv_cylinder)
        ev_list.extend(bg_out_fv_top_bottom)
        ev_list.extend(sn_ev)
        return ev_list

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
 #   ev_gen_ana = analytic_generator(M=2.0)
    sk = super_kamiokande(None, ev_gen)
    times = spectra.get_times()
    #times = np.linspace(0.01, 200, 20000)
    for time1, time2 in zip(times[1:], times[:-1]):
        sk.get_events(time2, time1)
