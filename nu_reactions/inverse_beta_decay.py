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

    def dcs_nue(self, Enu, E_target) -> float:
        return 0.0


    def dcs_cos_nue(self, Enu, cos):
        return np.zeros((len(cos), len(Enu[:,0]), len(Enu[:,1])))


    def cs_nue(self, Enu):
        return np.zeros(Enu.shape)


    def dcs_anue(self, Enu, E_target):

        dcsdE = np.zeros(Enu.shape, dtype=np.float_)

        Eth = ((self.__mn + self.__me)**2 - self.__mp**2)/2.0e0/self.__mp


        s =  np.where(Enu > Eth, self.__mp**2 + 2*self.__mp*Enu, self.__mp**2 + 2*self.__mp*Eth)
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
        dcsdE = np.where((E_target > 0.0), dsigmadE0\
        *(1.0e0+self.__alpha/np.pi*(6.0e0+1.5e0*np.log(0.5e0*self.__mp/np.abs(E_target)) + 1.2e0*np.sqrt((self.__me/np.abs(E_target))**3)))\
        ,0.0)

        dcsdE = np.where(Enu > Eth,dcsdE, 0.0)
        return dcsdE


    def dcs_cos_anue(self, Enu, cos):
        Enu_3d = np.tile(Enu, (len(cos), 1, 1))
        cos_new = np.reshape(cos, (len(cos),1,1))

        e = Enu_3d/self.__mp
        k = (1.0e0 + e)**2 - (e*cos_new)**2
        d = (self.__mn**2 - self.__mp**2 - self.__me**2)/2.0e0/self.__mp

        root2 = (Enu_3d-d)**2 - self.__me**2*k
        root2_new = np.where(root2 > 0.0, root2, 0.0)

        root = np.sqrt(root2_new)

        Ee = ((Enu_3d-d)*(1.0e0+e)+e*cos_new*root)/k
        pe = np.sqrt(np.abs(Ee**2 - self.__me**2))

        dcs_dE = self.dcs_anue(Enu, Ee)

        dcs_cos = pe*e/(1.0e0 + e*(1.0e0 - Ee/pe*cos_new))*dcs_dE
        dcs_cos = np.where(root2 < 0.0, 0.0, dcs_cos)
        return dcs_cos


    def cs_anue(self, Enu):

        n=1000
        cos_list = np.linspace(-1.0,1.0,n)
        cos1 = cos_list[:-1]
        cos2 = cos_list[1:]
        dcs_cos1 = self.dcs_cos_anue(Enu, cos1)
        dcs_cos2 = self.dcs_cos_anue(Enu, cos2)

        cs = np.sum(0.5*(dcs_cos1 + dcs_cos2)*(cos2[:,None,None] - cos1[:,None,None]), axis=0)
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


    def dcs_cos_nux(self, Enu, cos) -> float:
        return np.zeros((len(cos), len(Enu[:,0]), len(Enu[:,1])))


    def cs_nux(self, Enu):
        return np.zeros(Enu.shape)
    

if __name__ == "__main__":
    ibd = inverse_beta_decay()
    print(ibd.cs_anue(20))
