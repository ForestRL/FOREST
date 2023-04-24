from nu_reaction import nu_reaction
import numpy as np

class inverse_beta_decay(nu_reaction):
    """ Inverse beta decay cross-section given by Strumia and Vissani"""
    def __init__(self):
        self.__me = 0.510998e0 # electron mass
        self.__mp = 938.272e0 # proton mass
        self.__mn = 939.565e0 # neutron mass
        self.__delta = 1.293e0 # mass difference of the proton and the neutron
        self.__GF = 1.16637e-5*1.0e-6 # GF
        self.__costhetac = 0.9746e0
        self.__M = 938.9e0
        self.__g1_0 = -1.270e0
        self.__Mv2 = 0.71e6
        self.__MA2 = 1.0e6
        self.__xi = 3.706e0
        self.__hbar = 6.58212e-22
        self.__c = 2.997924e10
        self.__alpha = 1.0e0 / 137.036e0
        self.__mpi = 139.0e0


    def get_reaction_name(self) -> str:
        return "ibd_vissani"

    def dcs_nue(self, Enu:float, E_target:float) -> float:
        return 0.0


    def dcs_cos_nue(self, Enu:float, cos:float) -> float:
        return 0.0


    def cs_nue(self, Enu:float) -> float:
        return 0.0

    def dcs_anue(self, Enu, E_target) -> float:
            
        Eth = ((self.__mn + self.__me)**2 - self.__mp**2)/2.0e0/self.__mp

        if(Enu < Eth):
            dsigmadE = 0.0e0
            return 0.0

        s = self.__mp**2 + 2*self.__mp*Enu
        u = -(-s + 2*self.__mp*(Enu + E_target) - self.__me**2)
        t = self.__mn**2 - self.__mp**2 - 2.0e0*self.__mp*(Enu - E_target)

        f1 = (1-(1+self.__xi)*t/4/(self.__M**2))/(1-t/4/(self.__M**2))/((1-t/(self.__Mv2))**2)
        f2 = self.__xi/(1-t/4/(self.__M**2))/((1-t/(self.__Mv2))**2)
        g1 = self.__g1_0/((1-t/self.__MA2)**2)
        g2 = 2.0e0*self.__M**2*g1/(self.__mpi**2-t)

        At = (t-self.__me**2) * (0.25e0 *f1**2*(4.0e0*self.__M**2 + t + self.__me**2)
            +0.25e0 *g1**2*(-4.0e0*self.__M**2 + t + self.__me**2)
            +6.25e-2*f2**2*(t**2/self.__M**2 + 4.0e0*(t+self.__me**2))
            +0.25e0 *g2**2*t*self.__me**2/self.__M**2
            +0.5e0*f1*f2*(2.0e0*t + self.__me**2)  +  g1*g2*self.__me**2)\
            - self.__delta**2 * ( 6.25e-2*(4.0e0*f1**2 + f2**2*t/self.__M**2)\
            *( 4.0e0*self.__M**2 + t - self.__me**2)\
            +0.25e0 *g1**2*( 4.0e0*self.__M**2 - t + self.__me**2)\
            +0.25e0 *g2**2*(t - self.__me**2)*self.__me**2/self.__M**2\
            +0.5e0*f1*f2*(2.0e0*t - self.__me**2)  +  g1*g2*self.__me**2)\
            - 2.0e0*self.__me**2*self.__M*self.__delta*g1*(f1+f2)

        Bt = t*g1*(f1+f2)\
            + 0.25e0*self.__me**2*self.__delta*(f2**2+f1*f2+2.0e0*g1*g2)/self.__M
    
        Ct = 0.25e0*(f1**2+g1**2) - 6.25e-2*f2**2*t/self.__M**2

        Matrix = At - (s - u)*Bt + (s - u)**2*Ct

        dsigmadt = self.__GF**2*self.__costhetac**2/(2*np.pi*(s - self.__mp**2)**2)*Matrix
        dsigmadE0 = 2*self.__mp*dsigmadt*(self.__hbar*self.__c)**2

        dsigmadE = dsigmadE0\
        *(1.0e0+self.__alpha/np.pi*(6.0e0+1.5e0*np.log(0.5e0*self.__mp/E_target) + 1.2e0*np.sqrt((self.__me/E_target)**3)))

        return dsigmadE


    def dcs_cos_anue(self, Enu, cos) -> float:
        e = Enu/self.__mp
        k = (1.0e0 + e)**2 - (e*cos)**2
        d = (self.__mn**2 - self.__mp**2 - self.__me**2)/2.0e0/self.__mp

        root2 = (Enu-d)**2 - self.__me**2*k
        if(root2 < 0.0):
            dcs_cos = 0.0
            return dcs_cos
        
        root = np.sqrt(root2)

        Ee = ((Enu-d)*(1.0e0+e)+e*cos*root)/k
        pe = np.sqrt(Ee**2 - self.__me**2)

        dcs_dE = self.dcs_anue(Enu, Ee)

        dcs_cos = pe*e/(1.0e0 + e*(1.0e0 - Ee/pe*cos))*dcs_dE
        return dcs_cos


    def cs_anue(self, Enu:float) -> float:
        cs = 0

        n=1000
        cos_list = np.linspace(-1.0,1.0,n)
        for cos1, cos2 in zip(cos_list[:-1], cos_list[1:]):
            dcs_cos1 = self.dcs_cos_anue(Enu, cos1)
            dcs_cos2 = self.dcs_cos_anue(Enu, cos2)        
            cs += 0.5*(dcs_cos1 + dcs_cos2)*(cos2 - cos1)

        return cs


    def cs_table_anue(self, n_cos:int=200, n_e:int=300, E_min:float=0.0, E_max:float=300) -> list[np.array]:
        csdE = np.empty(n_e)
        dcsdEdcos = np.empty(n_e, n_cos)

        energy_list = np.linspace(E_min, E_max, n_e)
        cos_list = np.linspace(-1.0, 1.0, n_cos)

        for i, energy in enumerate(energy_list):
            csdE[i] = self.cs_anue(energy)
            for j, cos in enumerate(cos_list):
                pass



    def dcs_nux(self, Enu:float, E_target:float) -> float:
        return 0.0


    def dcs_cos_nux(self, Enu:float, cos:float) -> float:
        return 0.0


    def cs_nux(self, Enu:float) -> float:
        return 0.0
    

if __name__ == "__main__":
    ibd = inverse_beta_decay()
    print(ibd.cs_anue(20))
