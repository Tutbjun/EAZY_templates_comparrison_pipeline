##################### REBIN FUNCTION #####################
import numpy as np
def rebin(*args):
    return np.interp(*args)
def rebin(x_og, y_og, x_new, y_err=None):
    """ 
    Rebin the spectrum to get higher SN
    """
    og_sort = np.argsort(x_og)
    x_og = x_og[og_sort]
    y_og = y_og[og_sort]
    if y_err is not None: y_err = y_err[og_sort]

    x_new_offset = np.array([x_new[0]-(x_new[1] - x_new[0])/2] + [x_new[i] + (x_new[i+1] - x_new[i])/2 for i in range(len(x_new)-1)])
    x_og_offset = np.array([x_og[0]-(x_og[1] - x_og[0])/2] + [x_og[i] + (x_og[i+1] - x_og[i])/2 for i in range(len(x_og)-1)])
    
    y_new = np.zeros(len(x_new), dtype=np.float64)
    y_err_new = np.zeros(len(x_new), dtype=np.float64)
    for i in range(len(x_new)):
        index_within = np.where((x_og >= x_new_offset[i]) & (x_og <= x_new_offset[min(i+1, len(x_new_offset)-1)]))
        if len(index_within[0]) == 0:
            #interpolate
            y_new[i] = np.interp(x_new[i], x_og, y_og)#!not correct
            if y_err is not None:
                y_err_new[i] = np.interp(x_new[i], x_og, y_err)
            continue
        #y_new[i] = np.nanmean(y_og[index_within])#! does not factor in the width of the bin
        y_new[i] = np.trapz(y_og[index_within], x_og[index_within]) / (x_og[index_within][-1] - x_og[index_within][0])
        if y_err is not None:
            #y_err_new[i] = np.sqrt(np.sum(y_err[index_within]**2)) / len(index_within[0])#! does not factor in the width of the bin
            bindwidths_old = np.array([x_og_offset[i+1] - x_og_offset[i] for i in index_within[0]])
            y_err_new[i] = np.sqrt(np.sum(y_err[index_within]**2 * bindwidths_old) / (x_og[index_within][-1] - x_og[index_within][0]))
    return y_new, y_err_new
    
##################### REBIN FUNCTION #####################
def rebin(x_og, y_og, x_new, y_err):
    """ 
    Rebin the spectrum to get higher SN
    """
    #xerothly sort by x
    og_sort = np.argsort(x_og)
    x_og = x_og[og_sort]
    y_og = y_og[og_sort]
    y_err = y_err[og_sort]
    #firstly oversample the spectrum
    oversampleRate = 100
    y_new = np.zeros(len(x_new), dtype=np.float64)
    y_err_new = np.zeros(len(x_new), dtype=np.float64)
    """plt.errorbar(x_og, y_og, yerr=y_err, fmt='o')
    plt.show()"""
    x_oversample = np.linspace(x_og[0]/2-x_og[1]/2, x_og[-1]*3/2-x_og[-2]/2, len(x_og)*oversampleRate)
    y_oversample = np.interp(x_oversample, x_og, y_og)
    y_err_oversample = np.interp(x_oversample, x_og, y_err)
    y_err_oversample_resp = y_err_oversample**-1
    """plt.scatter(x_oversample, y_oversample, s=1)
    plt.plot(x_oversample, y_oversample+y_err_oversample, color='r')
    plt.plot(x_oversample, y_oversample-y_err_oversample, color='r')
    plt.plot(x_oversample, y_oversample+y_err_oversample_resp, color='b')
    plt.plot(x_oversample, y_oversample-y_err_oversample_resp, color='b')
    plt.ylim(min(y_oversample)-0.1, max(y_oversample)+0.1)
    plt.show()"""
    #now rebin the spectrum
    x_new_offset = np.array([x_new[0]-(x_new[1] - x_new[0])/2] + [x_new[i] + (x_new[i+1] - x_new[i])/2 for i in range(len(x_new)-1)])
    for i_n in range(len(x_new)):
        #first find the rebinvalues
        index_within = [i for i in range(len(x_oversample)) if (x_oversample[i] >= x_new_offset[i_n]) & (x_oversample[i] <= x_new_offset[min(i_n+1, len(x_new_offset)-1)])]
        if len(index_within) == 0:
            #interpolate
            y_new[i_n] = np.interp(x_new[i_n], x_oversample, y_oversample)
            y_err_new[i_n] = np.interp(x_new[i_n], x_oversample, y_err_oversample)
            continue
        #now rebin
        y_within = y_oversample[index_within]
        y_err_within = y_err_oversample[index_within]
        y_err_resp_within = y_err_oversample_resp[index_within]
        x_delta_within = x_oversample[index_within][-1] - x_oversample[index_within][0]
        y_new[i_n] = np.sum(y_within * x_delta_within * y_err_resp_within**2) / np.sum(x_delta_within * y_err_resp_within**2)#!Please check this formula
        #and then for error
        #y_err_new[i_n] = np.sqrt(np.sum(y_err_within**2 * x_delta_within * y_err_resp_within) / np.sum(x_delta_within * y_err_resp_within))
        y_err_new[i_n] = np.sqrt(np.sum(y_err_within**2) / len(y_err_within))#!Please check this formula
    """plt.errorbar(x_new, y_new, yerr=y_err_new, fmt='o')
    plt.show()"""
        

    return y_new, y_err_new

#do mock rebin
"""x_og = np.linspace(0, 100, 100)
y_og = np.random.normal(0, 1, len(x_og))
y_err = np.random.uniform(0, 1, len(x_og))
x_new = np.linspace(0, 100, 10)
y_new, y_err_new = rebin(x_og, y_og, x_new, y_err)"""

"""def rebin_spec(x, y, ye, xlo=None, xhi=None, bw=1):
    from spectres import spectres
    x_new = np.arange(xlo+bw, xhi, bw)
    y_new, ye_new = spectres(x_new, x, y, spec_errs=ye)
    return x_new, y_new, ye_new"""

from numba import jit

@jit(nopython=True)
def bandpass(wave_spec,flux_spec,wave_filt,flux_filt):#!now with numpy implementation
    #flux_filt = flux_filt / np.sum(flux_filt)#normalize integral
    #convert all wave to AA
    #wave_spec = np.array([wave_spec[i].to(u.AA).value for i in range(len(wave_spec))])
    #wave_filt = np.array([wave_filt[i].to(u.AA).value for i in range(len(wave_filt))])
    #use numpy interp
    filterrespons_interp = np.interp(wave_spec,wave_filt,flux_filt)#,left=0,right=0)
    for i in range(len(filterrespons_interp)):
        if wave_filt[i] < wave_spec[0] or wave_filt[i] > wave_spec[-1]:
            filterrespons_interp[i] = 0
    if sum(filterrespons_interp) == 0: return np.nan#*u.uJy
    filterrespons_interp = filterrespons_interp / np.sum(filterrespons_interp)
    #calculate product
    product = filterrespons_interp * flux_spec
    #integrate
    #integral = np.trapz(product,wave_spec)#/abs(wave_spec[-1] - wave_spec[0])
    integral = np.sum(product)
    #integral = integral / abs(wave_filt[-1] - wave_filt[0])
    return integral