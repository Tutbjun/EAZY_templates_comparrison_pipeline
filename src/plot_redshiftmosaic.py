import yaml
import os
##################### IMPORT SETTINGS #####################
#from settings.yaml
with open(os.path.join(os.path.dirname(__file__),"settings.yaml"), 'r') as stream:
    settings = yaml.safe_load(stream)
for key in settings['CONFIG']:
    globals()[key] = settings['CONFIG'][key]
for key in settings['MOSAICZ']:
    globals()[key] = settings['MOSAICZ'][key]

#known keys
zmin_global = zmin
zMax = zMax
mosTiling = mosTiling
mm = eval(mm)
figwidth = eval(figwidth)
DPI = DPI
ro = ro
zCharacteristic = zCharacteristic
SEM = SEM
sigmas = sigmas




def rnd_p(x):#, dmod=None):
    #print(x)
    def rnd(x:float, significantDigit=None):
        #print("ROunDING")
        if significantDigit is None:
            significantDigit = np.floor(np.log10(np.abs(x))) - 1
        rounded = round(x, -int(significantDigit))#, int(significantDigit)
        #if rounded == int(rounded): return int(rounded)
        #to string
        rounded = str(rounded)
        if significantDigit < 0:
            #print(rounded.split("."))
            #print(rounded)
            #print([rounded.split(".")[1][:-int(significantDigit)]])
            rounded = ".".join([rounded.split(".")[0]] + [rounded.split(".")[1][:-int(significantDigit)]])
            #print(rounded)
            #print([rounded.split(".")])
            #print(rounded)
            leadingDigits = rounded.split(".")[1]
            while leadingDigits.startswith("0") and not leadingDigits.endswith("0"): leadingDigits = leadingDigits[1:]
            while len(leadingDigits) < -significantDigit and str(leadingDigits)!='':
    
                #print("loop " + leadingDigits)
                rounded += "0"
                leadingDigits = rounded.split(".")[1]
                while leadingDigits.startswith("0") and not len(np.unique(leadingDigits))==1: leadingDigits = leadingDigits[1:]
            #print(rounded)
        return rounded
    #if its a sci_val, use the error for rounding
    if type(x) == sci_val:
        try:
            float(x.pm)
            assert not np.isnan(x.pm)
            assert not np.isinf(x.pm)
            assert not np.isnan(x.x)
            assert not np.isinf(x.x)
            err = x.pm
            X = x.x
            n = np.floor(np.log10(np.abs(err)))
            #print(err, n)
            rounded = float(rnd(X, n))
            err_rounded = float(rnd(err, n))
            return sci_val(rounded, err_rounded)
        except:
            if np.isnan(x.x) or np.isinf(x.x):
                return x
            return rnd(x.x)
    try:
        float(x)
        int(x)
    except:
        try:
            float(x[0])
            int(x[0])
            r = [rnd(i) for i in x]
            return r#, s
        except:
            return x#, None
    #if dmod is not None:
    #    #n = np.floor(np.log10(np.abs(err)))
    #    return rnd(x, -n)
    return rnd(x)


#TODO: needs tiling color by accumulative flux

etasStartIndex = 0
etas = []

from copy import copy, deepcopy
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

def clean_phots(z_spec, z_phot, zmin=5):
    mask = (z_spec > zmin) & (~np.isnan(z_spec))
    return z_spec[mask], z_phot[mask], mask

def get_deltaZ(z_spec, z_phot):
    return (z_phot - z_spec)/(1 + z_spec)

def get_outliers(z_spec, z_phot):
    deltaZs = get_deltaZ(z_spec, z_phot)
    return (deltaZs > 0.15) | (deltaZs < -0.15)

def get_bias(z_spec, z_phot):
    return np.mean(get_deltaZ(z_spec, z_phot))

def get_scatter(z_spec, z_phot):
    return np.std(get_deltaZ(z_spec, z_phot))

def get_eta(z_spec, z_phot):
    return np.sum(get_outliers(z_spec, z_phot))/len(z_spec)

def get_error(z_spec, z_phot):#the proportion of z_phot that is not fitted (<0 or nan)
    mask_havezpec = ~np.isnan(z_spec) & (z_spec >= 0)
    z_phot = z_phot[mask_havezpec]
    mask_havezphot = ~np.isnan(z_phot) & (z_phot >= 0)
    return np.sum(~mask_havezphot)/len(z_phot)

def get_stats(z_spec, z_phot):
    if z_spec is None or z_phot is None:
        return {
            'deltaZ': np.nan,
            'outliers': np.nan,
            'bias': np.nan,
            'scatter': np.nan,
            'eta': np.nan,
            'error': np.nan
        }
    stats = {}
    stats['deltaZ'] = get_deltaZ(z_spec, z_phot)
    stats['outliers'] = get_outliers(z_spec, z_phot)
    stats['bias'] = get_bias(z_spec, z_phot)
    stats['scatter'] = get_scatter(z_spec, z_phot)
    stats['eta'] = get_eta(z_spec, z_phot)
    stats['error'] = get_error(z_spec, z_phot)
    return stats

def get_stats_organized(df_out, ftempl_strs, ftempl_labels, zmin=5):
    deltaZ_s = {'default':[], 'modified':[]}
    deltaZ_all_s = {'default':[], 'modified':[]}
    outliers_s = {'default':[], 'modified':[]}
    outliers_all_s = {'default':[], 'modified':[]}
    bias_s = {'default':[], 'modified':[]}
    bias_all_s = {'default':[], 'modified':[]}
    scatter_s = {'default':[], 'modified':[]}
    scatter_all_s = {'default':[], 'modified':[]}
    etas_s = {'default':[], 'modified':[]}
    etas_all_s = {'default':[], 'modified':[]}
    error_s = {'default':[], 'modified':[]}
    error_all_s = {'default':[], 'modified':[]}
    for mode in ['default', 'modified']:
        mode_ftempl_lbls = copy(ftempl_labels)
        mode_ftempl_strs = copy(ftempl_strs)
        if mode == 'modified': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" in s]
        if mode == 'default': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" not in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" not in s]
        if len(mode_ftempl_lbls) == 0: continue
        for i, ftempl in enumerate(mode_ftempl_strs):
            if ftempl not in df_out.keys():
                continue
            #redshiftTbl = df_out[ftempl]
            redshiftTbl = {
                'z_spec': np.array(df_out[ftempl]['z_spec']),
                'z_phot': np.array(df_out[ftempl]['z_phot'])
            }
            redChar = {
                'z_spec': np.array(redshiftTbl['z_spec']),
                'z_phot': np.array(redshiftTbl['z_phot']),
            }
            
            redChar['z_spec'], redChar['z_phot'], _ = clean_phots(redChar['z_spec'], redChar['z_phot'], zmin=zmin)
            redshiftTbl['z_spec'], redshiftTbl['z_phot'], _ = clean_phots(redshiftTbl['z_spec'], redshiftTbl['z_phot'], zmin=zmin_global)

            stats_char = get_stats(redChar['z_spec'], redChar['z_phot'])
            stats_all = get_stats(redshiftTbl['z_spec'], redshiftTbl['z_phot'])
            deltaZ_s[mode].append(np.array([sci_val(stats_char['deltaZ'][i], np.nan) for i in range(len(stats_char['deltaZ']))]))
            deltaZ_all_s[mode].append(np.array([sci_val(stats_all['deltaZ'][i], np.nan) for i in range(len(stats_all['deltaZ']))]))
            outliers_s[mode].append(stats_char['outliers'])
            outliers_all_s[mode].append(stats_all['outliers'])
            bias_s[mode].append(sci_val(stats_char['bias'], np.nan))
            bias_all_s[mode].append(sci_val(stats_all['bias'], np.nan))
            scatter_s[mode].append(sci_val(stats_char['scatter'], np.nan))
            scatter_all_s[mode].append(sci_val(stats_all['scatter'], np.nan))
            etas_s[mode].append(sci_val(stats_char['eta'], np.nan))
            etas_all_s[mode].append(sci_val(stats_all['eta'], np.nan))
            error_s[mode].append(sci_val(stats_char['error'], np.nan))
            error_all_s[mode].append(sci_val(stats_all['error'], np.nan))
    return deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, bias_all_s, scatter_s, scatter_all_s, etas_s, etas_all_s, error_s, error_all_s

from numba import jit

percentileVal = None
def percentileFormula(sigmas):
    global percentileVal
    if percentileVal is not None:
        return percentileVal
    #evaluate:
    #1/sqrt(2pi)*integrate(exp(-x^2/2),0,sigmas)
    import scipy.integrate as spi
    percentileVal = spi.quad(lambda x: np.exp(-x**2/2), -sigmas, sigmas)[0]/np.sqrt(2*np.pi)*100/2
    return percentileFormula(sigmas)

#@jit(nopython=True)
def bootstrap_SEM(data, N=10000):
    nan_mask = np.isnan(data)
    data = data[~nan_mask]
    pN = len(data)
    x_avg_bs = np.empty(N)
    random_indices_all = np.random.randint(pN, size = (N,pN))# np.random.randint(N, size = N) creates N random numbers from 0 to N-1, i.e. N random indices
    data_resamplings = np.empty((N,pN))
    data_resamplings = data[random_indices_all]
    #data_resamplings = np.empty((1,pN))
    #data_resamplings[0] = data
    x_avg_bs = np.mean(data_resamplings, axis=1)
    data_new = data_resamplings.flatten()
    median = np.median(x_avg_bs)
    #from sigma setting, get lower and upper percent bounds
    ptil = percentileFormula(sigmas)
    if SEM:
        sigma1 = np.percentile(x_avg_bs, 50-ptil)
        sigma2 = np.percentile(x_avg_bs, 50+ptil)
    else:
        sigma1 = np.percentile(data_new, 50-ptil)
        sigma2 = np.percentile(data_new, 50+ptil)
    #print(rf"median: ${median}^{{+{sigma2-median}}}_{{-{median-sigma1}}}$")
    if median > 0.5:
        plotHistogramsOfRepeats(x_avg_bs, 'Distribution of the means')
    return median, median-sigma1, sigma2-median

def plotHistogramsOfRepeats(data, title):
    plt.hist(data, bins=100)
    plt.title(title)
    plt.show()

def get_stats_organized_mean(df_out, ftempl_strs, ftempl_labels, zmin=5):
    df_out_merged = get_resampleExtended(df_out)
    return get_stats_organized(df_out_merged, ftempl_strs, ftempl_labels, zmin=zmin)#TODO: validate


def get_stats_organized_variation(df_out, ftempl_strs, ftempl_labels, zmin=5):#TODO:validate
    #df_out_merged = get_resampleExtended(df_out)
    """for key in df_out.keys():
        
        df_out[key] = get_resampleExtended(df_out[key])"""
    deltaZ_s = {'default':[], 'modified':[]}
    deltaZ_all_s = {'default':[], 'modified':[]}
    outliers_s = {'default':[], 'modified':[]}
    outliers_all_s = {'default':[], 'modified':[]}
    bias_s = {'default':[], 'modified':[]}
    bias_all_s = {'default':[], 'modified':[]}
    scatter_s = {'default':[], 'modified':[]}
    scatter_all_s = {'default':[], 'modified':[]}
    etas_s = {'default':[], 'modified':[]}
    etas_all_s = {'default':[], 'modified':[]}
    error_s = {'default':[], 'modified':[]}
    error_all_s = {'default':[], 'modified':[]}
    for mode in ['default', 'modified']:
        deltaZ_char_part = {}
        deltaZ_all_part = {}
        outliers_char_part = {}
        outliers_all_part = {}
        bias_char_part = {}
        bias_all_part = {}
        scatter_char_part = {}
        scatter_all_part = {}
        etas_char_part = {}
        etas_all_part = {}
        error_char_part = {}
        error_all_part = {}
        mode_ftempl_lbls = copy(ftempl_labels)
        mode_ftempl_strs = copy(ftempl_strs)
        if mode == 'modified': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" in s]
        if mode == 'default': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" not in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" not in s]
        if len(mode_ftempl_lbls) == 0: continue
        for ftempl in mode_ftempl_strs:
            deltaZ_char_part[ftempl] = []
            deltaZ_all_part[ftempl] = []
            outliers_char_part[ftempl] = []
            outliers_all_part[ftempl] = []
            bias_char_part[ftempl] = []
            bias_all_part[ftempl] = []
            scatter_char_part[ftempl] = []
            scatter_all_part[ftempl] = []
            etas_char_part[ftempl] = []
            etas_all_part[ftempl] = []
            error_char_part[ftempl] = []
            error_all_part[ftempl] = []
        for i in range(max([len(df_out[key]) for key in df_out.keys()])):
            df_out_part = {}
            for key in df_out.keys():
                if len(df_out[key]) > i:
                    df_out_part[key] = df_out[key][i]
                else:
                    df_out_part[key] = None
            
            mask_all = None
            mask_char = None
            def stats_mask_calc(mask_all=None, mask_char=None):#!needs to have blanks for templates without data
                masks_all = []
                masks_char = []
                stats = {'char':[], 'all':[]}
                for j, ftempl in enumerate(mode_ftempl_strs):
                    if df_out_part[ftempl] is None:
                        stats_char = get_stats(None, None)
                        stats_all = get_stats(None, None)
                        stats['char'].append(stats_char)
                        stats['all'].append(stats_all)
                        continue
                    redshiftTbl = {
                        'z_spec': np.array(df_out_part[ftempl]['z_spec']),
                        'z_phot': np.array(df_out_part[ftempl]['z_phot'])
                    }
                    if mask_all is None:
                        redshiftTbl['z_spec'], redshiftTbl['z_phot'], mask = clean_phots(redshiftTbl['z_spec'], redshiftTbl['z_phot'], zmin=zmin_global)
                        masks_all.append(mask)
                    else:
                        redshiftTbl['z_spec'], redshiftTbl['z_phot'] = redshiftTbl['z_spec'][mask_all], redshiftTbl['z_phot'][mask_all]
                    redChar = {
                        'z_spec': np.array(redshiftTbl['z_spec']),
                        'z_phot': np.array(redshiftTbl['z_phot']),
                    }
                    if mask_char is None:
                        redChar['z_spec'], redChar['z_phot'], mask = clean_phots(redChar['z_spec'], redChar['z_phot'], zmin=zmin)
                        masks_char.append(mask)
                    else:
                        redChar['z_spec'], redChar['z_phot'] = redChar['z_spec'][mask_char], redChar['z_phot'][mask_char]
                    stats_char = get_stats(redChar['z_spec'], redChar['z_phot'])
                    stats_all = get_stats(redshiftTbl['z_spec'], redshiftTbl['z_phot'])
                    stats['char'].append(stats_char)
                    stats['all'].append(stats_all)
                return stats, masks_all, masks_char
            _, masks_all, masks_char = stats_mask_calc(mask_all, mask_char)
            #and-opperate masks
            mask_all = np.all(masks_all, axis=0)
            mask_char = np.all(masks_char, axis=0)
            stats, _, _ = stats_mask_calc(mask_all, mask_char)
            for j in range(len(stats['char'])):
                stats_char = stats['char'][j]
                stats_all = stats['all'][j]
                ftempl = mode_ftempl_strs[j]
                nanResult = False
                for key in stats_char.keys():
                    all_nan = np.all(np.isnan(stats_char[key]))
                    if all_nan:
                        nanResult = True
                        break
                if not nanResult:
                    deltaZ_char_part[ftempl].append(stats_char['deltaZ'])
                    deltaZ_all_part[ftempl].append(stats_all['deltaZ'])
                    outliers_char_part[ftempl].append(stats_char['outliers'])
                    outliers_all_part[ftempl].append(stats_all['outliers'])
                    bias_char_part[ftempl].append(stats_char['bias'])
                    bias_all_part[ftempl].append(stats_all['bias'])
                    scatter_char_part[ftempl].append(stats_char['scatter'])
                    scatter_all_part[ftempl].append(stats_all['scatter'])
                    etas_char_part[ftempl].append(stats_char['eta'])
                    etas_all_part[ftempl].append(stats_all['eta'])
                    error_char_part[ftempl].append(stats_char['error'])
                    error_all_part[ftempl].append(stats_all['error'])

        for ftempl in mode_ftempl_strs:
            deltaZ_char_part[ftempl] = np.array(deltaZ_char_part[ftempl])
            deltaZ_all_part[ftempl] = np.array(deltaZ_all_part[ftempl])
            outliers_char_part[ftempl] = np.array(outliers_char_part[ftempl])
            outliers_all_part[ftempl] = np.array(outliers_all_part[ftempl])
            bias_char_part[ftempl] = np.array(bias_char_part[ftempl])
            bias_all_part[ftempl] = np.array(bias_all_part[ftempl])
            scatter_char_part[ftempl] = np.array(scatter_char_part[ftempl])
            scatter_all_part[ftempl] = np.array(scatter_all_part[ftempl])
            etas_char_part[ftempl] = np.array(etas_char_part[ftempl])
            etas_all_part[ftempl] = np.array(etas_all_part[ftempl])
            error_char_part[ftempl] = np.array(error_char_part[ftempl])
            error_all_part[ftempl] = np.array(error_all_part[ftempl])
        deltaZ_char_val = {}
        deltaZ_all_val = {}
        outliers_char_val = {}
        outliers_all_val = {}
        bias_char_val = {}
        bias_all_val = {}
        scatter_char_val = {}
        scatter_all_val = {}
        etas_char_val = {}
        etas_all_val = {}
        error_char_val = {}
        error_all_val = {}
        for ftempl in mode_ftempl_strs:
            deltaZ_temp = deltaZ_char_part[ftempl].T
            deltaZ_char_val[ftempl] = np.array([sci_val(np.mean(deltaZ_temp[i]), np.std(deltaZ_temp[i])) for i in range(len(deltaZ_temp))])
            deltaZ_temp = deltaZ_all_part[ftempl].T
            deltaZ_all_val[ftempl] = np.array([sci_val(np.mean(deltaZ_temp[i]), np.std(deltaZ_temp[i])) for i in range(len(deltaZ_temp))])
            #outliers_temp = outliers_char_part[ftempl].T#!how does one calculate an error on a mask??
            if len(outliers_char_part[ftempl]) > 0: outliers_char_val[ftempl] = outliers_char_part[ftempl][0]
            else: outliers_char_val[ftempl] = np.array([])
            #sci_val(np.mean(outliers_temp), np.std(outliers_temp))
            #outliers_temp = outliers_all_part[ftempl].T
            #outliers_all_val[ftempl] = sci_val(np.mean(outliers_temp), np.std(outliers_temp))
            if len(outliers_all_part[ftempl]) > 0: outliers_all_val[ftempl] = outliers_all_part[ftempl][0]
            else: outliers_all_val[ftempl] = np.array([])
            AM, sigma1, sigma2 = bootstrap_SEM(bias_char_part[ftempl])
            bias_char_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(bias_all_part[ftempl])
            bias_all_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(scatter_char_part[ftempl])
            scatter_char_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(scatter_all_part[ftempl])
            scatter_all_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(etas_char_part[ftempl])
            etas_char_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(etas_all_part[ftempl])
            etas_all_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(error_char_part[ftempl])
            error_char_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
            AM, sigma1, sigma2 = bootstrap_SEM(error_all_part[ftempl])
            error_all_val[ftempl] = sci_val(AM, sigma1=sigma1, sigma2=sigma2)
        for ftempl in mode_ftempl_strs:
            deltaZ_s[mode].append(deltaZ_char_val[ftempl])
            deltaZ_all_s[mode].append(deltaZ_all_val[ftempl])
            outliers_s[mode].append(outliers_char_val[ftempl])
            outliers_all_s[mode].append(outliers_all_val[ftempl])
            bias_s[mode].append(bias_char_val[ftempl])
            bias_all_s[mode].append(bias_all_val[ftempl])
            scatter_s[mode].append(scatter_char_val[ftempl])
            scatter_all_s[mode].append(scatter_all_val[ftempl])
            etas_s[mode].append(etas_char_val[ftempl])
            etas_all_s[mode].append(etas_all_val[ftempl])
            error_s[mode].append(error_char_val[ftempl])
            error_all_s[mode].append(error_all_val[ftempl])
    pass
    return deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, bias_all_s, scatter_s, scatter_all_s, etas_s, etas_all_s, error_s, error_all_s

    #return get_stats_organized(df_out, ftempl_strs, ftempl_labels, zmin=zmin)#TODO: implement



class sci_val():
    def __format__(self, fmt):
        if fmt == '':
            return str(self)
        else:
            return float(self).__format__(fmt)
    def __repr__(self):
        return str(self.r) + "±" + str(self.pmr)
    def __init__(self, x, pm=None, sigma1=None, sigma2=None, npm=1):
        assert npm >= 1
        self.x = x
        self.pm = pm
        self.r = None
        if pm is not None:
            self.__round_pm__(pm,x,npm)
        elif sigma1 is not None and sigma2 is not None:
            self.__round_Q__(x,sigma1,sigma2,npm)
        else:
            self.r = x
            self.pmr = pm
            self.n = npm
    def __float__(self):
        return self.x
    def __round_pm__(self,pm,x, npm):
        self.sigma1 = None
        self.sigma2 = None
        try:
            #pm = self.ṕm
            float(pm)
            assert not np.isnan(pm)
            n = np.floor(np.log10(np.abs(pm)))-npm+1
            self.n = -n
            self.r = round(x, -int(n))
            self.pmr = round(pm, -int(n))
        except:
            self.r = x
            self.pmr = pm
            self.n = npm
    def __round_Q__(self,x,sigma1,sigma2, npm):
        self.pm = None
        self.pmr = None
        try:
            float(sigma1)
            float(sigma2)
            assert not np.isnan(sigma1)
            assert not np.isnan(sigma2)
            n = np.floor(np.log10(np.abs(min(sigma1,sigma2))))-npm+1
            self.n = -n
            self.r = round(x, -int(n))
            self.sigma1 = round(sigma1, -int(n))
            self.sigma2 = round(sigma2, -int(n))
        except:
            self.r = x
            self.sigma1 = sigma1
            self.sigma2 = sigma2
            self.n = npm
            
    def __str__(self):
        if self.pm is not None:
            r_str = str(self.r)
            if len(r_str.split(".")[-1]) > self.n:
                r_str = float(r_str)
                if not np.isnan(r_str) and not np.isinf(r_str):
                    if self.n <= 0: r_str = int(r_str)
                    else: r_str = round(r_str, int(self.n))
                r_str = str(r_str)
                while len(r_str.split(".")[-1]) < self.n:
                    r_str += "0"
            pmr_str = str(self.pmr)
            if len(pmr_str.split(".")[-1]) > self.n:
                pmr_str = float(pmr_str)
                if not np.isnan(pmr_str) and not np.isinf(pmr_str):
                    if self.n <= 0: pmr_str = int(pmr_str)
                    else: pmr_str = round(pmr_str, int(self.n))
                pmr_str = str(pmr_str)
                while len(pmr_str.split(".")[-1]) < self.n:
                    pmr_str += "0"
            return r_str + "±" + pmr_str
        elif self.sigma1 is not None and self.sigma2 is not None:
            r_str = str(self.r)
            if len(r_str.split(".")[-1]) > self.n:
                r_str = float(r_str)
                if not np.isnan(r_str) and not np.isinf(r_str):
                    if self.n <= 0: r_str = int(r_str)
                    else: r_str = round(r_str, int(self.n))
                r_str = str(r_str)
                while len(r_str.split(".")[-1]) < self.n:
                    r_str += "0"
            sigma1_str = str(self.sigma1)
            if len(sigma1_str.split(".")[-1]) > self.n:
                sigma1_str = float(sigma1_str)
                if not np.isnan(sigma1_str) and not np.isinf(sigma1_str):
                    if self.n <= 0: sigma1_str = int(sigma1_str)
                    else: sigma1_str = round(sigma1_str, int(self.n))
                sigma1_str = str(sigma1_str)
                while len(sigma1_str.split(".")[-1]) < self.n:
                    sigma1_str += "0"
            sigma2_str = str(self.sigma2)
            if len(sigma2_str.split(".")[-1]) > self.n:
                sigma2_str = float(sigma2_str)
                if not np.isnan(sigma2_str) and not np.isinf(sigma2_str):
                    if self.n <= 0: sigma2_str = int(sigma2_str)
                    else: sigma2_str = round(sigma2_str, int(self.n))
                sigma2_str = str(sigma2_str)
                while len(sigma2_str.split(".")[-1]) < self.n:
                    sigma2_str += "0"
            return r_str + "^{+" + sigma2_str + "}_{-" + sigma1_str + "}"
    def procent(self):
        if self.pm is not None and ~np.isnan(self.pm):
            x = self.x*100
            r = self.r*100
            pm = self.pm*100
            pmr = self.pmr*100
            n = self.n -2
            tempself = sci_val(x, pm)
            tempself.n = n
            return "["  + str(tempself) + "]%"
        elif self.sigma1 is not None and self.sigma2 is not None:
            x = self.x*100
            r = self.r*100
            sigma1 = self.sigma1*100
            sigma2 = self.sigma2*100
            n = self.n -2
            tempself = sci_val(x, sigma1=sigma1, sigma2=sigma2)
            tempself.n = n
            return "["  + str(tempself) + "]%"

from tqdm import tqdm
import pickle

"""def opperateWithBuffer(function,buffername,args):#open pickle if present, else create it
    tempdir = os.listdir('temp')
    if buffername not in tempdir:
        with open(os.path.join('temp',buffername), 'wb') as f:
            pickle.dump(function(*args), f)
    return pickle.load(open(os.path.join('temp',buffername), 'rb'))"""

def get_resampleExtended(df_out_resamplings):#!TODO HERE: dont take a bootstrapped average, just concatinate them all!
    new_df_out = {}
    for key in df_out_resamplings.keys():
        #avegare z_phot and z_spec
        z_phots = []
        z_specs = []
        new_df_out[key] = {}
        for i in range(len(df_out_resamplings[key])):
            z_phots.append(df_out_resamplings[key][i]['z_phot'])
            z_specs.append(df_out_resamplings[key][i]['z_spec'])
        z_phots = np.array(z_phots)
        z_specs = np.array(z_specs)
        #replace all negative values with nan
        for i in range(len(z_phots)):
            negative_mask = (z_phots[i] < 0) | (z_specs[i] < 0) | (np.isnan(z_phots[i])) | (np.isnan(z_specs[i]))
            z_phots[i][negative_mask] = np.nan
            z_specs[i][negative_mask] = np.nan
            #z_specs = z_specs.T
            #z_phots = z_phots.T
        #use bootstrap to get the mean and error
        
        """if f'boostrapbuffer_{key}.pkl' not in os.listdir('temp'):
            z_phot_means = []
            z_spec_means = []
            for i in tqdm(range(len(z_phots)), desc="Calculating bootstrap averages for " + key, total=len(z_phots)):
                val, _, _ = bootstrap_SEM(z_phots[i])
                z_phot_means.append(val)
                val, _, _ = bootstrap_SEM(z_specs[i])
                z_spec_means.append(val)
            with open(os.path.join('temp',f'boostrapbuffer_{key}.pkl'), 'wb') as f:
                pickle.dump((z_phot_means, z_spec_means), f)
        else:
            z_phot_means, z_spec_means = pickle.load(open(os.path.join('temp',f'boostrapbuffer_{key}.pkl'), 'rb'))
        z_phot_means = np.array(z_phot_means).T
        z_spec_means = np.array(z_spec_means).T
        new_df_out[key]['z_phot'] = z_phot_means
        new_df_out[key]['z_spec'] = z_spec_means"""
        z_specs = z_specs.flatten()#!new
        z_phots = z_phots.flatten()
        new_df_out[key]['z_phot'] = z_phots
        new_df_out[key]['z_spec'] = z_specs
    return new_df_out

from matplotlib.colors import LogNorm

def plot(df_out,df_out_resamplings,ftempl_strs,ftempl_labels,runTime):
    global etasStartIndex
    global etas
    etasStartIndex = 0
    etas = []
    #precalc stats
    deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, _, scatter_s, _, etas_s, _, _, _, have_full_resample = get_data(df_out,df_out_resamplings,ftempl_strs,ftempl_labels,runTime)

    for mode in ['default', 'modified']:
        mode_ftempl_lbls = copy(ftempl_labels)
        mode_ftempl_strs = copy(ftempl_strs)
        if mode == 'modified': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" in s]
        if mode == 'default': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" not in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" not in s]
        if len(mode_ftempl_lbls) == 0: continue
        
        mosaicLen = ceil(len(mode_ftempl_lbls)/mosTiling)
        
        figsize = (1*figwidth,1*figwidth/mosTiling*mosaicLen)
        tiling = (mosaicLen,round(mosTiling))
        if len(mode_ftempl_lbls) < mosTiling:
            figsize = (1*figwidth/mosTiling*len(mode_ftempl_lbls),1*figwidth/mosTiling)
            tiling = (1,len(mode_ftempl_lbls))

        figMos, axsMos = plt.subplots(*tiling, figsize=figsize, dpi=DPI, facecolor=(1,1,1,0))
        #transpose and flatten
        if len(mode_ftempl_lbls) > 1:
            #axsMos = axsMos.T
            axsMos = axsMos.flatten()
            axsBlank = axsMos[len(mode_ftempl_lbls):]
            axsMos = axsMos[:len(mode_ftempl_lbls)]
        else:
            axsMos = [axsMos]
            axsBlank = []
        plt.subplots_adjust(wspace=0, hspace=0)
        if len(axsBlank) > 0:
            for ax in axsBlank:
                ax.axis('off')

        #color the backgrounds of the plots with colormap RdYlGn by eta
        temp = [etas_s[k] for k in etas_s.keys()]
        etas = []
        for i in range(len(temp)): etas += temp[i]
        etas = np.array([e.x for e in etas])
        
        maxEr, minEr = np.min(1-etas),np.max(1-etas)#np.min(etas)**(-ro), np.max(etas)**(-ro)
        norm = mpl.colors.Normalize(vmin=maxEr, vmax=minEr)
        cmap = mpl.cm.RdYlGn
        cmap.set_under(color=cmap(maxEr))
        for i in range(etasStartIndex,len(axsMos)+etasStartIndex):
            j = i - etasStartIndex
            color = cmap(norm((1-etas[i])))#**(-ro)))
            color = (*color[:-1], 0.1)#set alpha = 0.1
            axsMos[j].set_facecolor(color)
        etasStartIndex = len(axsMos)
        print(f"Have {len(df_out[mode_ftempl_strs[0]]['z_spec'])} objects")
        if have_full_resample:
            print(f"Have resamples, gathering full z range")
            df_out_2use = get_resampleExtended(df_out_resamplings)
        else:
            df_out_2use = df_out
        for i, ftempl in enumerate(mode_ftempl_strs):
            if ftempl not in df_out.keys():
                continue
            
            redshiftTbl = df_out_2use[ftempl]
            redChar = {
                'z_spec': np.array(redshiftTbl['z_spec']),
                'z_phot': np.array(redshiftTbl['z_phot']),
            }
            #mask = (redChar['z_spec'] > zCharacteristic) & (~np.isnan(redChar['z_spec'])) & (~np.isnan(redChar['z_phot']) & (redChar['z_phot'] > 0))
            #redChar['z_spec'] = redChar['z_spec'][mask]
            #redChar['z_phot'] = redChar['z_phot'][mask]
            redChar['z_spec'], redChar['z_phot'], _ = clean_phots(redChar['z_spec'], redChar['z_phot'], zmin=zCharacteristic)
            
            deltaZ = deltaZ_s[mode][i]
            deltaZ_all = deltaZ_all_s[mode][i]
            outliers = outliers_s[mode][i]
            outliers_all = outliers_all_s[mode][i]
            bias = bias_s[mode][i]
            scatter = scatter_s[mode][i]
            eta = etas_s[mode][i]
            
            #Remove sub zero and nan
            x = np.array(redshiftTbl['z_spec'])
            y = np.array(redshiftTbl['z_phot'])
            
            
            x, y, _ = clean_phots(x, y, zmin=0)

            #calculate chi2 with linear regression of x and y
            chi2 = np.sum((y - x)**2/1)/len(x)#/len(x)#!dunno what to put for sigma

            xmin, xmax = np.c_[[x, y]].min(), np.c_[[x, y]].max()

            underCharc = (x < zCharacteristic)
            x_under, y_under = x[underCharc], y[underCharc]
            #hexbins with log
            axsMos[i].hexbin(x_under, y_under, gridsize=70, mincnt=1, edgecolors='none', cmap='viridis', alpha=0.2, extent=[xmin, xmax, xmin, xmax],zorder=200, norm=LogNorm())
            axsMos[i].hexbin(redChar['z_spec'], redChar['z_phot'], gridsize=70, mincnt=1, edgecolors='none', cmap='viridis', extent=[xmin, xmax, xmin, xmax], norm=LogNorm())

            #plt.colorbar(axsMos[i].collections[0], ax=axsMos[i], label='log$_{10}$(N)')
            axsMos[i].plot([0, 100], [0, 100], c='r', ls='--', lw=0.5)
            axsMos[i].plot([zCharacteristic, zCharacteristic], [0, 14.5], c='k', ls='--', lw=0.3)

            #plot the 15% and 85% quantiles
            zeropoint = [0,0]
            mindevL = [0,-0.15]
            mindevH = [0,0.15]
            maxdevL = [100,(-0.15-100*0.15+100)]
            maxdevH = [100,(0.15+100*0.15+100)]
            axsMos[i].plot([mindevL[0],maxdevL[0]], [mindevL[1],maxdevL[1]], c='r', ls='--', lw=0.5)
            axsMos[i].plot([mindevH[0],maxdevH[0]], [mindevH[1],maxdevH[1]], c='r', ls='--', lw=0.5)
            print(mindevL, maxdevL)
            print(mindevH, maxdevH)
            #okay just plot something within frame to show that it works
            #axsMos[i].plot([1, 1], [5,4], c='r', ls='-', lw=0.5)
            
            #dict_stat = hmod.phot_spec_zs_stats(y, x)

            #annotate in top left
            annotAnchor = (0.05, 0.95)
            #axsMos[i].set_title(f'{ftempl}', fontsize=10)
            ax, ay = 50, 95        
            axsMos[i].annotate(f'{mode_ftempl_lbls[i]}', xy=annotAnchor, xycoords='axes fraction', fontsize=10, ha='left', va='top', xytext=(ax-48, ay), textcoords='axes points')

            #stats annotate
            sx, sy = 45, -55.5
            etastr = (eta.procent() + "$").replace("%", "\%")
            axsMos[i].annotate(f'$\eta={etastr}', xy=annotAnchor, xycoords='axes fraction', fontsize=5.8, ha='right', va='top', xytext=(ax+sx, ay+sy), textcoords='axes points')
            #etas = np.append(etas, eta**(-1))
            if mode_ftempl_lbls[i] == 'F45k' or mode_ftempl_lbls[i] == 'F60k': asterix = '*'
            else: asterix = ''
            #axsMos[i].annotate(f'$\eta_{{0+}}$={round(outliers_all.sum()/len(outliers_all),2):.2f}{asterix}', xy=annotAnchor, xycoords='axes fraction', fontsize=8, ha='left', va='top', xytext=(ax+sx, ay+sy-10), textcoords='axes points')
            vals = np.array([v.x for v in deltaZ[~outliers]])
            errs = np.array([v.pm for v in deltaZ[~outliers]])
            #sigma = np.median(np.abs(vals - np.median(vals)))*1.4826
            #sigma_err = np.median(np.abs(errs - np.median(errs)))*1.4826
            #sigma = sci_val(sigma, np.nan)
            #sigmastr = sigma.procent().replace("nan", "?")
            sigmastr = (scatter.procent() + "$").replace("%", "\%")
            biasstr = (bias.procent() + "$").replace("%", "\%")
            axsMos[i].annotate(f'$\sigma_{{nmad}}={sigmastr}', xy=annotAnchor, xycoords='axes fraction', fontsize=5.8, ha='right', va='top', xytext=(ax+sx, ay+sy-10), textcoords='axes points')
            axsMos[i].annotate(f'$\Delta_{{bias}}={biasstr}', xy=annotAnchor, xycoords='axes fraction', fontsize=5.8, ha='right', va='top', xytext=(ax+sx, ay+sy-20), textcoords='axes points')


            axsMos[i].set_xlim(zmin_global,zMax)
            axsMos[i].set_ylim(0,zMax)

            axsMos[i].xaxis.set_minor_locator(MultipleLocator(1))
            axsMos[i].yaxis.set_minor_locator(MultipleLocator(1))
            axsMos[i].xaxis.set_major_locator(MultipleLocator(4))
            axsMos[i].yaxis.set_major_locator(MultipleLocator(4))

            if i % mosTiling != 0: axsMos[i].set_yticklabels([])
            else: axsMos[i].set_ylabel('$z_{phot}$', fontsize=10)
            if i < len(mode_ftempl_lbls)-mosTiling: axsMos[i].set_xticklabels([])
            else: axsMos[i].set_xlabel('$z_{spec}$', fontsize=10)

        for i in range(len(axsMos)):
            if i > len(mode_ftempl_lbls)-1:
                axsMos[i].axis('off')

        
        

        """if len(ftempl_strs) % 3 != 0:
            axsMos[-1].axis('off')
        if len(ftempl_strs) % 3 == 1:
            axsMos[-2].axis('off')"""

        #axis label
        #figMos.text(0.5, 0.05, '$z_{spec}$', ha='center', va='center')
        #figMos.text(0.05, 0.5, '$z_{phot}$', ha='center', va='center', rotation='vertical')
        #plt.tight_layout()
        #plt.show()
        figMos.savefig(f'./figures/forpaper/zs_mosaic_{mode}_{runTime}.png', dpi=DPI, bbox_inches='tight', facecolor=(1,1,1,0))
        figMos.clf()
        plt.close(figMos)
        #show by loading image
        img = plt.imread(f'./figures/forpaper/zs_mosaic_{mode}_{runTime}.png')
        fig = plt.figure(figsize=(1*figwidth,1*figwidth/mosTiling*mosaicLen), dpi=200, facecolor=(1,1,1,0))
        ax = fig.add_axes([0,0,1,1])
        ax.axis('off')
        ax.imshow(img)

import pandas as pd

def get_data(df_out,df_out_resamplings,ftempl_strs,ftempl_labels,runTime):
    deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, bias_all_s, scatter_s, scatter_all_s, etas_s, etas_all_s, error_s, error_all_s = get_stats_organized(df_out,ftempl_strs,ftempl_labels, zmin=zCharacteristic)
    deltaZ_s_res_m, deltaZ_all_s_res_m, outliers_s_res_m, outliers_all_s_res_m, bias_s_res_m, bias_all_s_res_m, scatter_s_res_m, scatter_all_s_res_m, etas_s_res_m, etas_all_s_res_m, error_s_res_m, error_all_s_res_m = get_stats_organized_mean(df_out_resamplings,ftempl_strs,ftempl_labels, zmin=zCharacteristic)
    _, _, _, _, bias_s_res_var, bias_all_s_res_var, scatter_s_res_var, scatter_all_s_res_var, etas_s_res_var, etas_all_s_res_var, error_s_res_var, errors_all_s_res_var = get_stats_organized_variation(df_out_resamplings,ftempl_strs,ftempl_labels, zmin=zCharacteristic)
    have_full_resample = True
    for key in df_out_resamplings.keys():
        if len(df_out_resamplings[key]) == 0:
            have_full_resample = False
            break
    for mode in ['default', 'modified']:
        for i in range(len(deltaZ_all_s_res_m[mode])):
            if len(deltaZ_all_s_res_m[mode][i]) > 0:
                deltaZ_s[mode][i] = deltaZ_s_res_m[mode][i]
                deltaZ_all_s[mode][i] = deltaZ_all_s_res_m[mode][i]
                outliers_s[mode][i] = outliers_s_res_m[mode][i]
                outliers_all_s[mode][i] = outliers_all_s_res_m[mode][i]
                bias_s[mode][i] = bias_s_res_m[mode][i]
                bias_all_s[mode][i] = bias_all_s_res_m[mode][i]
                scatter_s[mode][i] = scatter_s_res_m[mode][i]
                scatter_all_s[mode][i] = scatter_all_s_res_m[mode][i]
                etas_s[mode][i] = etas_s_res_m[mode][i]
                etas_all_s[mode][i] = etas_all_s_res_m[mode][i]
                error_s[mode][i] = error_s_res_m[mode][i]
                error_all_s[mode][i] = error_all_s_res_m[mode][i]
        for j in range(len(bias_s[mode])):
            bias_s[mode][j] = sci_val(bias_s[mode][j].x, sigma1=bias_s_res_var[mode][j].sigma1, sigma2=bias_s_res_var[mode][j].sigma2)
        
            bias_all_s[mode][j] = sci_val(bias_all_s[mode][j].x, sigma1=bias_all_s_res_var[mode][j].sigma1, sigma2=bias_all_s_res_var[mode][j].sigma2)
        
            scatter_s[mode][j] = sci_val(scatter_s[mode][j].x, sigma1=scatter_s_res_var[mode][j].sigma1, sigma2=scatter_s_res_var[mode][j].sigma2)
            scatter_all_s[mode][j] = sci_val(scatter_all_s[mode][j].x, sigma1=scatter_all_s_res_var[mode][j].sigma1, sigma2=scatter_all_s_res_var[mode][j].sigma2)
        
            etas_s[mode][j] = sci_val(etas_s[mode][j].x, sigma1=etas_s_res_var[mode][j].sigma1, sigma2=etas_s_res_var[mode][j].sigma2)
            etas_all_s[mode][j] = sci_val(etas_all_s[mode][j].x, sigma1=etas_all_s_res_var[mode][j].sigma1, sigma2=etas_all_s_res_var[mode][j].sigma2)
        
            error_s[mode][j] = sci_val(error_s[mode][j].x, sigma1=error_s_res_var[mode][j].sigma1, sigma2=error_s_res_var[mode][j].sigma2)
            error_all_s[mode][j] = sci_val(error_all_s[mode][j].x, sigma1=errors_all_s_res_var[mode][j].sigma1, sigma2=errors_all_s_res_var[mode][j].sigma2)
    return deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, bias_all_s, scatter_s, scatter_all_s, etas_s, etas_all_s, error_s, error_all_s, have_full_resample

def table(df_out,df_out_resamplings,ftempl_strs,ftempl_labels,runTime):
    global etasStartIndex
    global etas
    etasStartIndex = 0
    etas = []
    #precalc stats
    deltaZ_s, deltaZ_all_s, outliers_s, outliers_all_s, bias_s, bias_all_s, scatter_s, scatter_all_s, etas_s, etas_all_s, errors_s, errors_all_s, have_full_resample = get_data(df_out,df_out_resamplings,ftempl_strs,ftempl_labels,runTime)
    for mode in ['default', 'modified']:
        mode_ftempl_lbls = copy(ftempl_labels)
        mode_ftempl_strs = copy(ftempl_strs)
        if mode == 'modified': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" in s]
        if mode == 'default': 
            mode_ftempl_lbls = [s for s in mode_ftempl_lbls if ";" not in s]
            mode_ftempl_strs = [s for s in mode_ftempl_strs if ";" not in s]
        if len(mode_ftempl_lbls) == 0: continue
        df = {}
        
            
        insertions = {
            #'deltaZ': 'deltaZ_s[{mode}]',
            #'deltaZ_all': 'deltaZ_all_s[{mode}]',
            #'outliers': 'outliers_s[{mode}]',
            #'outliers_all': 'outliers_all_s[{mode}]',
            '$\Delta_{bias,z>5}$': 'bias_s[{mode}]',
            '$\Delta_{bias,z>0.1}$': 'bias_all_s[{mode}]',
            '$\sigma_{nmad,z>5}$': 'scatter_s[{mode}]',
            '$\sigma_{nmad,z>0.1}$': 'scatter_all_s[{mode}]',
            '$\eta_{z>5}$': 'etas_s[{mode}]',
            '$\eta_{z>0.1}$': 'etas_all_s[{mode}]',
            '$e_{z>5}$': 'errors_s[{mode}]',
            '$e_{z>0.1}$': 'errors_all_s[{mode}]',
            #chi2
        }
        """for key in k
            df[key] = []"""
    #for mode in ['default', 'modified']:
        df = {}
        df['templ.'] = [
            f'${s}$' for s in mode_ftempl_lbls
        ]
        
        for key in insertions.keys():
            #x, digit = rnd_p(eval(insertions[key].format(mode=repr(mode))))
            df[key] = eval(insertions[key].format(mode=repr(mode)))
        print(df)

        df = pd.DataFrame(df)
        for key in df.keys():
            for i in range(len(df[key])):
                #df[key][i] = rnd_p(df[key][i])
                if type(df[key][i]) == sci_val:
                    """df[key][i].x *= 100
                    df[key][i].r *= 100
                    df[key][i].pm *= 100
                    df[key][i].pmr *= 100
                    df[key][i].n -= 2
                    df[key][i] = "[" + str(df[key][i]) + "]\%"#!potentielt lav øvre nedre +- formattering her"""
                    df[key][i] = df[key][i].procent().replace("%","\%")
                    df[key][i] = df[key][i].replace("[","$[")
                    df[key][i] = df[key][i].replace("%","%$")
        cmap = mpl.cm.RdYlGn
        colortable = []
        for key in list(df.keys())[1:]:
            if "±" in df[key][0]:
                values = np.array([float(s.split("±")[0].split("[")[-1]) for s in df[key]])
            elif "^" in df[key][0]:
                values = np.array([float(s.split("^")[0].split("[")[-1]) for s in df[key]])
            #errs = np.array([float(s.pm) for s in df[key]])
            min = np.min(np.abs(values.astype(float))**-1)
            max = np.max(np.abs(values.astype(float))**-1)
            colortable.append([])
            for i in range(len(df[key])):
                colortable[-1].append(np.array(cmap((np.abs(float(values[i]))**-1 - min)/(max-min))[:-1]))
                #colortable[-1][-1][-1] = 0.1
                colortable[-1][-1] *= 1.2
                colortable[-1][-1] = np.clip(colortable[-1][-1], 0, 1)
        colortable = np.array(colortable).transpose(1,0,2)

        #print(df)
        #save latex to txt file
        with open(f'./figures/forpaper/zs_table_temp_{runTime}.tex', 'w') as f:
            f.write(df.to_latex(index=False))
        texlines = []
        with open(f'./figures/forpaper/zs_table_temp_{runTime}.tex', 'r') as f:
            texlines = f.readlines()
        try:
            os.remove(f'./figures/forpaper/zs_table_temp_{runTime}.tex')
        except FileNotFoundError:
            pass
        for i in range(len(texlines)):
            if texlines[i].startswith("\\begin{tabular}"):
                texlines[i] = "\\begin{deluxetable*}{c|cccccccc}[t]\n"
            if texlines[i].startswith("\\toprule"):
                texlines[i] = "\\label{tab:zs_stats}\n\\tabletypesize{\\scriptsize}\n\\tablewidth{0pt}\n\\tablecaption{Photometric redshift performance metrics}\n"
                header = texlines[i+1]
                header = header.split(" \\\\\n")[0].split(" & ")
                texlines[i+1] = "\\tablehead{\n" + "&".join(["\\colhead{" + h + "}" for h in header]) + "} \\\n\\startdata\n"
            if texlines[i].startswith("\\midrule"):
                texlines[i] = ""
            if texlines[i].startswith("\\bottomrule"):
                texlines[i] = "\\enddata\n\\hline\n\\tablecomments{autogenerated...}\n\\end{deluxetable*}"
            if texlines[i].startswith("\\end{tabular}"):
                texlines[i] = ""
        dataReached = False
        print("COLORTABLE")
        print(colortable)
        di = 0
        
        for i in range(len(texlines)):
            #if i == 0 or i == 1 or i == 2: continue
            print(texlines[i])
            if not dataReached:
                if "\\startdata" in texlines[i]:
                    dataReached = True
                continue
            else:
                print("colouring line " + str(i))
                print(texlines[i])
                newLine = texlines[i].split(" & ")[0] + " & "
                dj = 0
                for j,val in enumerate(texlines[i].split(" & ")[1:]):
                    try:
                        print("c: " + str(colortable[di][dj]))
                        newLine += "\\cellcolor[HTML]{" + str(mpl.colors.rgb2hex(colortable[di][dj])).upper().replace("#","") + "} " + val + " & "
                    except IndexError:
                        newLine += val + " & "
                    dj += 1
                print(newLine)
                di += 1
                #texlines[i] = newLine + "\\\\\n"


        try:
            os.remove(f'./figures/forpaper/zs_table_{mode}_{runTime}.tex')
        except FileNotFoundError:
            pass
        with open(f'./figures/forpaper/zs_table_{mode}_{runTime}.tex', 'w') as f:
            f.writelines(texlines)
            f.close()