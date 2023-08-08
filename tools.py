import numpy as np

def get_value_random(x_list, pdf)->float:
    """
        returns a interpolated value of x_list at random with PDF.
    """
    
    if(len(x_list) != len(pdf)):
        return None

    if(np.sum(pdf) < 1.0e-100):
        return 0.0

    cum_y = np.cumsum(pdf)
    max_y = cum_y[-1]

    normal_cum_y = cum_y/max_y
    rand = np.random.rand()
#    print(rand, normal_cum_y)
    if(rand < normal_cum_y.min()):
        return x_list[0]

    i = np.argmax(normal_cum_y[normal_cum_y - rand < 0])
    value = (x_list[i+1] - x_list[i])/(normal_cum_y[i+1] - normal_cum_y[i])*(rand - normal_cum_y[i]) + x_list[i]
    return value 

def get_value_random_hist(low_edges, high_edges, histogram)->float:
    cum_hist = np.cumsum(histogram)
    cdf = cum_hist/cum_hist[-1]

    rand = np.random.rand()
    i = np.searchsorted(cdf, rand)
    low_edge = low_edges[i]
    high_edge = high_edges[i]
    value = (high_edge - low_edge)/(cdf[i] - cdf[i-1])*(rand - cdf[i-1]) + low_edge
    return value


def interpolate_dcsdcos_ene(ene, ene_bins_all, dcsdEdcos):
    """
        returns dcsdEdcos array with interplation at ene
    """

    i = np.argmax(ene_bins_all[ene_bins_all - ene < 0])
    new_dcsdcos = (dcsdEdcos[i+1,:] -dcsdEdcos[i,:])/(ene_bins_all[i+1] - ene_bins_all[i])\
                            *(ene - ene_bins_all[i]) + dcsdEdcos[i,:]
    return new_dcsdcos

def interpolate_1d_array(x, x_array, array_1d):
    """
        returns scalar from 1d array with interplation at time
    """

    if(x <= x_array.min()):
        return 0.0
    if(x >= x_array.max()):
        return 0.0
    
    i = np.argmax(x_array[x_array - x < 0])
    new_array = (array_1d[i+1] - array_1d[i])/(x_array[i+1] - x_array[i])\
                    *(x - x_array[i]) + array_1d[i]
    
    return new_array    

def interpolate_2d_array(x, x_array, array_2d):
    """
        returns 1d array from 2d array with interplation at time
    """

    if(x <= x_array.min()):
        return np.zeros_like(array_2d[0,:])
    if(x >= x_array.max()):
        return np.zeros_like(array_2d[-1,:])
    
    i = np.argmax(x_array[x_array - x < 0])
    new_array = (array_2d[i+1,:] - array_2d[i,:])/(x_array[i+1] - x_array[i])\
                    *(x - x_array[i]) + array_2d[i,:]
    return new_array    

def interpolate_3d_array(x, x_array, array_3d):
    """
        returns 2d array from 3d array with interplation at time
    """

    i = np.argmax(x_array[x_array - x < 0])
    new_array = (array_3d[i+1,:,:] - array_3d[i,:,:])/(x_array[i+1] - x_array[i])\
                    *(x - x_array[i]) + array_3d[i,:,:]
    return new_array


def fermi(kT:float, n:float, x_array):
    """
        returns the fermi function of order n
    """
    f = x_array**n/(1.+np.exp(x_array/kT))
    return f
