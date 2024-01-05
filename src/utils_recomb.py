#!/usr/bin/env python
# coding: utf-8

# In[45]:


#load spectres of templates at last point
import os
from astropy.table import Table
import numpy as np

import eazy
import matplotlib.pyplot as plt
from numba import jit
from astropy.io import fits
import helper_module as hmod
from astropy.table import Table, join
import eazy_routines as ez
from astropy import units as u
import time
from astropy.cosmology import Planck18
import pandas as pd
import shutil
import sys, os
import warnings
from scipy.optimize import curve_fit


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

useSquareInLoss = True

"""runTime = int(time.time())
__file__ = os.getcwd()+'/'+"gen-template-matrix_byZ.ipynb"
cosmo = Planck18
set2Limit = "SPHINX"""
#set2Limit = "EMlines"
#set2Limit = "fsps_hot/45k"


#templCnt = 4
#templCnt = 4
#templCnt = 39
#cCap = 0.5
"""spectrasetCap = 400
trainsetCap = 2000
trainsetFrac = 0.75#0.8
useSquareInLoss = True
matrixName = input("Optimized matrix name: ")
out_path = "templates-custom"#!temp
prevMatrixName = 'SPHINX_4_400c_0_279'#'EMLines_08p_3_14'#'serverCCap_0_28'
onlyTestset = False"""

#path = "templates-custom/45k"
#path = "templates-custom/EMlines"#!temp
#path = "templates-custom/SPHINX"
#out_path = "templates-custom/45k"

"""path = os.path.dirname(os.path.realpath(__file__)) + "/" + path
out_path = os.path.dirname(os.path.realpath(__file__)) + "/" + out_path

templates_in = [f for f in os.listdir(path) if "bin1.fits" in f or "bin0.fits" in f or ".spec" in f or ".fits" in f]
spectras = []
for temp in templates_in:
    if ".fits" in temp:
        tab = Table.read(path+"/"+temp)
        flux = np.asarray(tab["flux"]).astype(float).T#.T[0]
        wave = np.asarray(tab["wave"]).astype(float).T#.T[0]

    elif ".spec" in temp:
        tab = np.loadtxt(path+"/"+temp)
        flux = tab.T[-1]
        wave = tab.T[0]
    spectras.append(flux)
spectras = np.array(spectras)
for i in range(len(spectras)):
    spectras[i] = spectras[i] / np.sum(spectras[i])"""

"""if len(spectras) > spectrasetCap:
    pick = np.random.choice(range(len(spectras)), spectrasetCap, replace=False)
    spectras = spectras[pick]
    print("spectra capped")"""


#cCaps = np.full((templCnt, spectras.shape[0]), 1)

#flt = eazy.filters.FilterFile()

"""filts_nircam = {
        'F090W': 363,
        'F115W': 364,
        'F150W': 365,
        'F182M': 370,
        'F200W': 366,
        'F210M': 371,
        'F277W': 375,
        'F335M': 381,
        'F356W': 376,
        'F410M': 383,
        'F430M': 384,
        'F444W': 358,
        'F460M': 385,
        'F480M': 386
}

filts_HST = {
        'F105W': 202,
        'F125W': 203,
        'F140W': 204,
        'F160W': 205,
        'F435W': 233,
        'F606W': 214,
        'F775W': 216,
        'F814W': 239,
        'F850LP': 240
}"""

"""filts = {**filts_nircam, **filts_HST}

mw_reddening = ez.get_atten_dict(filts)

zps = [1.0]*len(filts)"""


"""inname = "hlsp_jades_jwst_nircam_goods-s-deep_photometry_v1.0_catalog_large_withSpec.fits"
inpath = os.path.join(os.getenv('astrodata'), 'gds/jades/phot', inname)"""

"""with fits.open(inpath) as hdul:
    print(hdul.info())"""

"""tab = Table.read(inpath, hdu=6)
tab_redshifts = Table.read(inpath, hdu=9)"""

"""ext = '_CIRC1'
cols_fluxes = hmod.get_matches(ext, tab.columns, exclude='_ei')
cols_f = np.sort(hmod.get_matches(ext, cols_fluxes, exclude='_e'))
cols_fe = np.sort(hmod.get_matches('_e', cols_fluxes))
cols_fluxes = list(np.vstack([cols_f, cols_fe]).T.flatten())
cols = list(np.insert(cols_fluxes, 0, ['ID', 'RA', 'DEC', 'z_spec','z_spec_quality']))

tab = join(tab, tab_redshifts['ID', 'z_spec', 'z_spec_quality'], join_type='inner', keys='ID')
tab_out = tab[cols]"""

"""keys = np.array(list(mw_reddening.keys()))
for c in cols_fluxes:
    tab_out[c].unit = u.nJy
    tab_out[c] = tab_out[c].to(u.uJy)
    
    # apply MW reddening
    matches = hmod.get_matches(keys, c, get_idxs=True)
    key = keys[np.int32(matches[:,0])][0]
    tab_out[c] *= mw_reddening[key]"""

#for c in cols_f:
#    cnew = c.replace(ext, '_flux')
#    tab_out.rename_column(c, cnew)

#for c in cols_fe:
#    cnew = c.replace(ext+'_e', '_err')
#    tab_out.rename_column(c, cnew)


#tab_out_Aqual = tab_out[tab_out['z_spec_quality'] == 'A']

#tab_out_Aqual_train = tab_out_Aqual[np.random.choice(range(len(tab_out_Aqual)), min(int(len(tab_out_Aqual)*trainsetFrac),trainsetCap), replace=False)]
#tab_out_Aqual_test = tab_out_Aqual[~np.isin(tab_out_Aqual['ID'], tab_out_Aqual_train['ID'])]

# save EAZY table
#os.makedirs('data', exist_ok=True)
#tab_out.write('./data/gds_jades_eazy.fits', format='fits', overwrite=True)
#tab_out_Aqual.write('./data/gds_jades_eazy_Aqual.fits', format='fits', overwrite=True)
#tab_out_Aqual_train.write('./data/gds_jades_eazy_Aqual_train.fits', format='fits', overwrite=True)
#tab_out_Aqual_test.write('./data/gds_jades_eazy_Aqual_test.fits', format='fits', overwrite=True)
def writeDummyTemplateset(size:int):#just writes a dummy version
    if "optimized_template" not in os.listdir("temp"):
        os.makedirs("temp/optimized_template")
        print("Created intermediate folder: "+"temp/optimized_template")
    
    with open(f'temp/optimized_template/optimized.param', 'w') as f:
        for i in range(size):
            f.write(f'1 temp/optimized_template/optimized_{i}.spec\n')
        f.close()
    
    if f"optimized_0.spec" not in os.listdir("temp/optimized_template"):
        X = np.linspace(91,1e8,10000)
        Y = np.ones(len(X))
        writeList = np.array([X, Y]).T
        for i in range(size):
            np.savetxt(f'temp/optimized_template/optimized_{i}.spec', writeList)
        


def writeTemplateset(ks,spectras,wave):
    newSpectras = performMatrixOperation(ks,spectras)
    #save the new templateset
    if "optimized_template" not in os.listdir("temp"):
        os.makedirs("temp/optimized_template")
        print("Created intermediate folder: "+"temp/optimized_template")
    with open(f'temp/optimized_template/optimized.param', 'w') as f:
        #plt.clf()
        for i,flux,_ in zip(range(len(newSpectras)),newSpectras,wave):
            if f"optimized_{i}.spec" in os.listdir("temp/optimized_template"):
                os.remove(f'temp/optimized_template/optimized_{i}.spec')
                #print("Removed intermediate file: "+f'{intermediateFolder}/{i}.spec')
            writeList = np.array([wave, flux]).T
            np.savetxt(f'temp/optimized_template/optimized_{i}.spec', writeList)
            plt.plot(wave, flux)
            #print("Created template file: "+f'{intermediateFolder}/{i}.spec')
            f.write(f'{i+1} temp/optimized_template/optimized_{i}.spec\n')
        f.close()
        #log axis
        #plt.xscale("log")
        #plt.yscale("log")
        #plt.xlim(100,100000)
        #plt.savefig("specs.png")
    print("updated template set")

"""intermediateFolder = "spectras_temp"
if not os.path.exists(intermediateFolder):
    os.makedirs(intermediateFolder)
    print("Created intermediate folder: "+intermediateFolder)
if os.path.exists(f'{intermediateFolder}/paramFile.param'):
    os.remove(f'{intermediateFolder}/paramFile.param')
    print("Removed intermediate file: "+f'{intermediateFolder}/paramFile.param')
with open(f'{intermediateFolder}/paramFile.param', 'w') as p:
    saveTemplateset(spectras[:templCnt])
    for i in range(templCnt):
        p.write(f'{i+1} {intermediateFolder}/{i}.spec\n')
p.close()"""


#cat_name_train = 'gds_jades_eazy_Aqual_train'
#cat_name_test = 'gds_jades_eazy_Aqual_test'
#cat_path_train = f'./data/{cat_name_train}.fits'
#cat_path_test = f'./data/{cat_name_test}.fits'
#keys_id = ['ID id', 'RA ra', 'DEC dec', 'z_spec z_spec']

#templ_paths = f'{intermediateFolder}/paramFile.param'



@jit(nopython=True)
def performMatrixOperation(matrix, spectra):
    return np.dot(matrix,spectra)

import utils_astro
def runEAZY(photZ, trainOrTest="train"):
    """if doTestset:
        cat_name = cat_name_test
        cat_path = cat_path_test
    else:
        cat_name = cat_name_train
        cat_path = cat_path_train
    tpath = f'{intermediateFolder}/paramFile.param'
    with open(tpath) as f:
        lines = f.readlines()
        f.close()
    opath = f'./eazy-output/{runTime}'
    params = {"cat_path": cat_path,
        "templ_path": tpath,
        "out_path": opath,
        "FIX_ZSPEC": 'n',
        "USE_ZSPEC_FOR_REST": 'n',
        "Z_MAX": 15,
        "H0": cosmo.H0,
        "OMEGA_M": cosmo.Om0,
        "OMEGA_L": cosmo.Ode0,
        "CATALOG_FORMAT": 'fits'}
    # write eazy config files
    ___, fnames = ez.write_config(cat_name, filts, zps, keys_id,
        out_path=opath)
    # run eazy"""
    print("################ EAZY START ###################")
    idx = None
    #from astropy.cosmology import Planck18 as cosmo
    """tpath = f'temp/optimized_template/optimized.param'
    opath = f'eazy-output/temporrrrrrrrrrrrrraaaarytest'
    cat_path = f'temp/train.fits'
    cat_name = 'train'
    filts = utils_astro.get_all_filters()
    zps = [1.0]*len(filts)
    keys_id = ['ID id', 'RA ra', 'DEC dec', 'z_spec z_spec']
    
    params = {"cat_path": cat_path,
        "templ_path": tpath,
        "out_path": opath,
        "FIX_ZSPEC": 'n',
        "USE_ZSPEC_FOR_REST": 'n',
        "Z_MAX": 15,
        "H0": cosmo.H0,
        "OMEGA_M": cosmo.Om0,
        "OMEGA_L": cosmo.Ode0,
        "CATALOG_FORMAT": 'fits'}
    # write eazy config files
    ___, fnames = ez.write_config(cat_name, filts, zps, keys_id,
        out_path=opath)"""
    """with HiddenPrints():
        #blockPrint()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")"""
            #_, __ = ez.run_eazy(params, fnames, n_proc=-1, idx=idx)
    #photZ_new = eazy.photoz.PhotoZ(params=photZ.params, translate_file=photZ.translate_file)
    #photZ.params['CATALOG_FILE'] = f'temp/{trainOrTest}.fits'
    
    photZ.iterate_zp_templates(n_proc=-1, verbose=False, idx=idx, update_zeropoints=False)#!dunno if its good that it doesnt update zeropoints
    photZ.fit_catalog(n_proc=-1, verbose=False, idx=idx)
    plt.clf()

    #photZ_new.fit_catalog(n_proc=-1, verbose=True, idx=idx)
    #photZ = photZ_new
    """zout, hdu = photZ.standard_output(idx=idx)
    print(zout)
            #warnings.resetwarnings()
        #enablePrint()
    #zout, hdu = 
    #read the output
    oname = 'eazy'
    cat_out_name = 'gds_jwst_nircam_STSCI'
    fpath = f'temp/eazy-output-optimizer/{trainOrTest}.zout.fits'
    tbl = Table.read(fpath)
    names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]
    df = tbl[names].to_pandas()
    #change df 'id' to 'ID'
    df.rename(columns={'id': 'ID'}, inplace=True)"""

    #merge with catalog zspec #!instead, pull in a catalogue with zspec and just merge? either that or it is already there?
    """fpath = './data/gds_jades_eazy_Aqual_train.fits'
    tbl = Table.read(fpath, hdu=1)    
    names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]
    names = [name for name in names if 
        'ID' in name or
        'EAZY_z_a' in name or 
        'z_spec_source' in name or 
        'z_spec_quality' in name or
        'z_spec_reference' in name
    ]
    df_spec = tbl[names].to_pandas()
    
    # add spec info to the photo table
    df = pd.merge(df, df_spec, on='ID', how='left')"""

    print("################ EAZY DONE ###################")
    return photZ

def meassureBadnessOfZ(matrix,spectras,waves,dset,trainOrTest="train"):
    #function that encapsulates the whole framework of testing EAZY on a set of spectras
    #templateSpectras = performMatrixOperation(matrix,spectras)
    writeTemplateset(matrix,spectras,waves)
    Zs = runEAZY(dset,trainOrTest=trainOrTest)
    """if "old_zphot.csv" in os.listdir():
        z_phot_old = np.loadtxt(f'old_zphot.csv', delimiter=',')
    else:
        z_phot_old = Zs.zbest"""
    z_spec = Zs.ZSPEC
    z_phot = Zs.zbest
    """np.savetxt(f'old_zphot.csv', z_phot, delimiter=',')
    if len(z_phot) == len(z_phot_old):
        print("change in z_phot: "+str(np.sum(np.abs(z_phot_old - z_phot))))
    else:
        print("eh...")"""
    zDeltas = z_spec - z_phot
    global useSquareInLoss
    if useSquareInLoss: return np.sum(zDeltas**2)/len(zDeltas), Zs
    else: return np.sum(np.abs(zDeltas))/len(zDeltas), Zs



#matrix = np.zeros((templCnt, spectras.shape[0]), dtype=np.float64)
#for i in range(templCnt):
#    matrix[i][i] = 1

#if prevMatrixName != '':
#    matrix = np.loadtxt(f'{out_path}/{prevMatrixName}.csv', delimiter=',')
#    print("Loaded previous matrix")

"""@jit(nopython=True)
def matrixOp(matrix, spectras):
    vec = np.zeros((templCnt, 1))
    for i in range(templCnt):
        vec[i] = np.sum(matrix[i] * spectras[i])
    return vec"""


#c = 0
#polynomium = lambda x, a, b, c: a*x**2 + b*x + c

def loadspecs(templatepaths:list):
    specs = []#! uhm, the resulting flux's are 14 long, wich is a thing that stems way back. why is this??
    waves = []
    for path in templatepaths:
        print(path)
        if ".fits" in path:
            tab = Table.read(path)
            flux = np.asarray(tab["flux"]).astype(float).T#.T[0]
            wave = np.asarray(tab["wave"]).astype(float).T#.T[0]

        elif ".spec" in path:
            tab = np.loadtxt(path)
            flux = tab.T[-1]
            wave = tab.T[0]
        specs.append(np.array(flux))
        waves.append(np.array(wave))
    
    return waves, specs


import pickle
import utils_math
def GETSPEC(template_path):
    if 'specs.pkl' in os.listdir('temp'):
        with open('temp/specs.pkl', 'rb') as f:
            specs = pickle.load(f)
            f.close()
        with open('temp/waves.pkl', 'rb') as f:
            waves = pickle.load(f)
            f.close()
        return waves[0], specs
    tmpl_folder = template_path.split('/')[0]
    spec_paths = np.loadtxt(template_path, dtype=str, usecols=1)
    #spec_paths = [os.path.join(tmpl_folder, "/".join(s.split('/')[1:])) for s in spec_paths]
    waves, specs = loadspecs(spec_paths)
    for i in range(len(specs)):
        if len(specs[i].shape) > 1:#!hot fix, we should figure out why there are so many sub-spectras in the first place
            specs[i] = specs[i][0]
    #if there are different waves, rebin the higher res ones to the lowest
    minWaveRes_index = np.argmin([len(w) for w in waves])
    minWave = waves[minWaveRes_index]
    minWaveRes = len(minWave)
    for i in range(len(waves)):
        if len(waves[i]) > minWaveRes:
            flux_new, _ = utils_math.rebin(waves[i], specs[i], minWave, y_err=np.ones(len(specs[i])), progress=True)
            specs[i] = flux_new
            waves[i] = minWave
    specs = np.array(specs)
    waves = np.array(waves)
    for i in range(len(specs)):
        if len(specs[i].shape) > 1:#!this normalization is super important to get right
            specs[i] = specs[i] / np.trapz(np.mean(specs[i],axis=0), waves[i])
        else:
            specs[i] /= np.trapz(specs[i], waves[i])
    #save specs
    with open('temp/specs.pkl', 'wb') as f:
        pickle.dump(specs, f)
        f.close()
    with open('temp/waves.pkl', 'wb') as f:
        pickle.dump(waves, f)
        f.close()
    return waves[0], specs

def INITSTATE(wave,specs, newSize, photZ_test):
    #check for existing state
    temp_folder = os.listdir("temp")
    if "recomb_state.pkl" in temp_folder:
        with open("temp/recomb_state.pkl", "rb") as f:
            state = pickle.load(f)
            assert len(state["cMaxs"]) >= 1
            assert len(state["ks"]) >= 1
            assert np.shape(state["cMaxs"]) == np.shape(state["ks"])
            assert len(state["cMaxs"]) == newSize
            assert len(state["cMaxs"][0]) == len(specs)
            assert len(state["losses_test"]) >= 1
            assert not np.isnan(state["ks"]).any()
            f.close()
    else:#genereate starting state
        ks_temp = np.diag(np.ones(newSize))
        #resample matrix by "stretching" horizontally
        ks = np.zeros((newSize, len(specs)))
        for i in range(len(specs)):
            ks[:,i] = ks_temp[int(i/len(specs)*newSize)]
        for i in range(len(ks)):
            ks[i] = ks[i] / np.sum(ks[i])
        #genereate cMaxs
        cMaxs = np.ones(ks.shape)
        state = {"cMaxs":cMaxs, "ks":ks, "losses_test":[(run:=meassureBadnessOfZ(ks,specs,wave,dset=photZ_test,trainOrTest="test"))[0]]}
        with open("temp/recomb_state.pkl", "wb") as f:
            pickle.dump(state, f)
            f.close()
    return state

from copy import deepcopy
def testC(c,N,matrix,spectras,wave,photZ_test):
    matrix_test = deepcopy(matrix + c*N)
    #plt.imshow(matrix_test, cmap="hot")
    #plt.savefig("matrix_test.png")
    return (run:=meassureBadnessOfZ(matrix_test,spectras,wave,dset=photZ_test,trainOrTest="test"))[0]

def backup_state(state):
    old_states = os.listdir("temp/oldstates")
    old_states = [s for s in old_states if "recomb_state" in s]
    old_states = np.sort(old_states)[::-1]
    for s in old_states:
        #move old states to backup folder
        oldness = int(s.split("_")[2].split(".")[0])
        new_oldness = oldness + 1
        os.rename(f"temp/oldstates/{s}", f"temp/oldstates/recomb_state_{new_oldness}.pkl")
    #move newest state to oldstates
    if "recomb_state.pkl" in os.listdir("temp"):
        os.rename("temp/recomb_state.pkl", "temp/oldstates/recomb_state_1.pkl")
    #save current state
    with open("temp/recomb_state.pkl", "wb") as f:
        pickle.dump(state, f)
        f.close()

def normalizeMatrix(matrix):
    for j in range(len(matrix)):
        for k in range(len(matrix[j])):
            if np.isnan(matrix[j][k]): matrix[j][k] = 0
            if np.isinf(matrix[j][k]): matrix[j][k] = 1
            if matrix[j][k] > 1: matrix[j][k] = 1
            if matrix[j][k] < 0: matrix[j][k] = 0
        if np.sum(matrix[j]) >= 1: 
            matrix[j] = matrix[j] / np.sum(matrix[j])
    return matrix

def FOUND_PROCESS(loss,losses_train,matrix,i,nCol,c,N):
    #global betterFound
    #global losses_train
    betterFound = True
    
    matrix = matrix + c*N
    matrix = normalizeMatrix(matrix)
    print("IMPROVEMENT FOUND!!!")
    np.savetxt(f'temp/optimized_template/matrix/matrix__{i}_{nCol}.csv', matrix, delimiter=',')
    print(rf"$log_{{10}}\Delta$ loss: {np.log10(np.abs(loss-losses_train[-1]))}")
    losses_train.append(loss)
    
    #losses_train.append(loss)
    return betterFound, losses_train, matrix

def skip_state(i,nCol,specCnt,inner):
    if inner: return skip_state_inner(i,nCol,specCnt)
    doContinue = False
    if os.path.exists(f'temp/_state.txt'):
        with open(f'temp/_state.txt', 'r') as f:
            state = f.read()
            f.close()
        state = state.split(" ")
        i_pred = int(state[0])
        nCol_pred = int(state[1])
        nCol_pred = nCol_pred
        if nCol_pred >= specCnt:
            nCol_pred = 0
            i_pred = i_pred + 1
        if i != i_pred:
            doContinue = True
        #if nCol != nCol_pred: continue
    return doContinue

def skip_state_inner(i,nCol,specCnt):
    doContinue = False
    if os.path.exists(f'temp/_state.txt'):
        with open(f'temp/_state.txt', 'r') as f:
            state = f.read()
            f.close()
        state = state.split(" ")
        i_pred = int(state[0])
        nCol_pred = int(state[1])
        nCol_pred = nCol_pred
        if nCol_pred < 0:
            nCol_pred = specCnt-1
            i_pred = i_pred + 1
        if i != i_pred: doContinue = True
        if nCol != nCol_pred: doContinue = True
    return doContinue

def write_state(i,nCol,specCnt):
    with open(f'temp/_state.txt', 'w') as f:
        if nCol <= 0:
            f.write(f"{i+1} {specCnt-1}")
        else:
            f.write(f"{i} {nCol-1}")
        f.close()

def eval_rnd_Npoint(cCaps, nCol, pickBuffer):
    p1,p2 = 0,0
    doContinue = False
    pickStack = np.argwhere(cCaps[:,nCol] == np.amax(cCaps[:,nCol])).flatten()
    foundPoint = False
    #print("eval_rnd_Npoint")
    while not foundPoint and len(pickStack) > 0:
        #print(f"pickStack: {pickStack}")
        np.random.seed(int(time.time()))
        randInt = np.random.randint(0,len(pickStack))
        #print(f"randInt: {randInt}")
        #p1 = np.random.choice(pickStack.flatten())
        p1 = pickStack[randInt]
        p2 = nCol
        #print(f"p1: {p1}, p2: {p2}")
        if [p1,p2] in pickBuffer:
            pickStack = np.delete(pickStack, np.argwhere(pickStack == p1))
            continue
        else: foundPoint = True
    if not foundPoint: 
        doContinue = True
        return doContinue, p1, p2
    return doContinue, p1, p2

def recombine_templates(photZ_train, photZ_test, template_path, newSize, filts, zps, keys_id, cosmo, maxZ=20, zmin=0, itters=100):
    #load template spectras
    wave, specs = GETSPEC(template_path)
    #load/make state of matrix and cMaxs
    #TODO: update photZ's with variable template set path
    state = INITSTATE(wave, specs, newSize, photZ_test)
    print("################### START #######################")
    c=0
    
    losses_train = []
    losses_test = []
    writeTemplateset(state["ks"],specs,wave)
    backup_state(state)
    losses_train.append((run:=meassureBadnessOfZ(state["ks"],specs,wave,dset=photZ_train,trainOrTest="train"))[0])
    
    try:
        for i in range(itters):
            if skip_state(i,0,len(specs),inner=False): continue
            print(f"################### loop i={i} #######################")
            betterFound = False
            pickBuffer = []
            if i == 0: 
                for j in range(len(specs)): 
                    for k in range(newSize):
                        pickBuffer.append([j,k])
            for nCol,NCnt in list(enumerate(range(len(specs))))[::-1]:#nesting of c-optimization
                #read state file
                if skip_state(i,nCol,len(specs),inner=True): continue
                #write state file
                write_state(i,nCol,len(specs))
                #read cCaps
                #if os.path.exists(f'{out_path}/cCaps.csv'):
                cCaps = state["cMaxs"]
                #write matrix
                print(f"################### loop nCol={nCol} #######################")
                #set N to a matrix with a random index =1
                betterFound = False
                N = np.zeros(state['ks'].shape,dtype=np.float64)
                cont, p1, p2 = eval_rnd_Npoint(cCaps, nCol, pickBuffer)
                if cont: continue
                N[p1][p2] = 1
                cCap = cCaps[p1][p2]
                cMin = (-state['ks'][p1][p2])*cCap
                cMax = (1-state['ks'][p1][p2])*cCap
                pickBuffer.append([p1,p2])
                plt.imshow(state['ks'], cmap="hot")
                plt.colorbar()
                plt.savefig(f"matrix_{i}.png")
                plt.savefig("matrix_current.png")
                plt.close()
                
                def updatePlt(lossByCFunc, poly=None):
                    figure, axs = plt.subplots(1, 1, figsize=(5,5))
                    lines = [
                        axs.plot([], [])[0],# polynomialfit
                        axs.plot([], [], marker=".", linestyle="")[0]# lossByCFunc
                    ]
                    axs.set_xlabel("c")
                    axs.set_ylabel("loss")
                    lines[1].set_data(np.array(lossByCFunc)[:,0], np.array(lossByCFunc)[:,1])
                    axs.set_ylim([losses_train[-1], (max(np.array(lossByCFunc)[:,1])-min(np.array(lossByCFunc)[:,1]))*1.1+min(np.array(lossByCFunc)[:,1])])
                    axs.set_xlim([min(np.array(lossByCFunc)[:,0])*1.1, max(np.array(lossByCFunc)[:,0])*1.1])
                    if poly is not None:
                        lines[0].set_data(np.linspace(-cCap,cCap,100), np.polyval(poly, np.linspace(-cCap,cCap,100)))
                    #figure.canvas.draw()
                    #plt.savefig(f"lossByCFunc_{i}_{NCnt}.png")
                    plt.savefig("lossByCFunc_current.png")
                    plt.close()

                lossByCFunc = [[0,losses_train[-1]]]
                updatePlt(lossByCFunc)
                ####################

                #initially test at pm cCap
                
                if cMax == 0:
                    c = cMin/2
                    print("cMin/2")
                else:
                    c = cMax
                    print("cMax")
                print(f"test @ cCap: {c}...")
                loss = testC(c,N,state['ks'],specs,wave,photZ_train)
                print("first c test done")
                if loss < losses_train[-1]:
                    cCnt = 1
                    #betterFound = True
                    #losses_train.append(loss)
                    #matrix = matrix + c*N
                    #matrix = normalizeMatrix(matrix)
                    betterFound, losses_train, state['ks'] = FOUND_PROCESS(loss,losses_train,state['ks'],i,nCol,c,N)
                    backup_state(state)

                    continue
                lossByCFunc.append([c,loss])
                updatePlt(lossByCFunc)
                if cMin == 0:
                    c = cMax/2
                    print("cMax/2")
                else:
                    c = cMin
                    print("cMin")
                print(f"test @ -cCap: {c}...")
                loss = testC(c,N,state['ks'],specs,wave,photZ_train)
                if loss < losses_train[-1]: 
                    cCnt = 2
                    #betterFound = True
                    #losses_train.append(loss)
                    #matrix = matrix + c*N
                    #matrix = normalizeMatrix(matrix)
                    betterFound, losses_train, state['ks'] = FOUND_PROCESS(loss,losses_train,state['ks'],i,nCol,c,N)
                    backup_state(state)
                    
                    continue
                lossByCFunc.append([c,loss])
                updatePlt(lossByCFunc)
                print("test at -cCap done")
                #find c's by solving polynomials for ))st likely bottom
                for cCnt in range(3,6):#loop of finding c's
                    print(f"################### loop cCnt={cCnt} #######################")
                    order = cCnt - 2
                    if order < 3: order = 2
                    #do polynomial fit
                    #poly = curve_fit(polynomium, np.array(lossByCFunc)[:,0], np.array(lossByCFunc)[:,1], p0=[1,0,1])[0]
                    #try:#TODO: traceback the nan/inf issue in lossByCFunc
                    poly = np.polyfit(np.array(lossByCFunc)[:,0], np.array(lossByCFunc)[:,1], deg=order)
                    print(f"poly: {poly}")
                    #insert dots and polyfit
                    updatePlt(lossByCFunc, poly)
                    #differentiate
                    polyDif = np.polyder(poly)
                    #find roots
                    roots = np.real(np.roots(polyDif))#!exclude complex roots
                    print(roots)
                    #pick a root at random but weighted by how isolated the root is and how low the loss is at that point
                    c = roots[np.argmin(np.polyval(poly, roots))]
                    rootVal = np.polyval(poly, c)
                    print(rootVal)
                    if c < cMin or c > cMax or rootVal > losses_train[-1]:#i <= 5 and c < 0: 
                        prob = (np.max(np.polyval(poly, np.linspace(0, cCap, 100)))-np.polyval(poly, np.linspace(0, cCap, 100)))
                        c = np.random.choice(np.linspace(0, cCap, 100), 1, p=prob/np.sum(prob))[0]
                    #except:
                    #    print("Plonomium exception!!!")
                    #    c = np.random.uniform(cMin,cMax)
                    """if abs(c) >= cCap:
                        prob = (np.max(np.polyval(poly, np.linspace(-cCap, cCap, 100)))-np.polyval(poly, np.linspace(-cCap, cCap, 100)))
                        c = np.random.choice(np.linspace(-cCap, cCap, 100), 1, p=prob/np.sum(prob))"""
                    
                    print(f"trying c: {c}")
                    loss = testC(c,N,state['ks'],specs,wave,photZ_train)
                    lossByCFunc.append([c,loss])
                    if loss < losses_train[-1]: 
                        #found(loss)
                        #betterFound = True
                        #matrix = matrix + c*N
                        #matrix = normalizeMatrix(matrix)
                        #losses_train.append(loss)
                        #np.savetxt(f'{out_path}/{matrixName}_{i}_{nCol}.csv', matrix, delimiter=',')
                        betterFound, losses_train, state['ks'] = FOUND_PROCESS(loss,losses_train,state['ks'],i,nCol,c,N)
                        backup_state(state)
                        break
                    plt.close()
                    plt.imshow(state['ks'] + c*N, cmap="hot")
                    plt.colorbar()
                    plt.savefig(f"matrix_{i}.png")
                    plt.savefig("matrix_current.png")
                    plt.close()
                    plt.imshow(N, cmap="hot")
                    plt.colorbar()
                    plt.savefig(f"N_{i}.png")
                    print(f"loss: {loss}, c: {c}")
                updatePlt(lossByCFunc,poly)
                np.savetxt(f'temp/optimized_template/matrix/matrix__{i}_{nCol}.csv', state['ks'], delimiter=',')
                #"""if betterFound: 
                if not betterFound: 
                    cCaps[p1][p2] = cCap/2
                if betterFound: 
                    losses_train.append(loss)
                    #if nCol == 0:
                    #    losses_test.append((run:=meassureBadnessOfZ(state['ks'],specs,wave,dset=photZ_test,trainOrTest="test"))[0])
                    cCaps[p1][p2] = cCap
                #write cCaps
                state["cMaxs"] = cCaps
                
                #state["ks"] = matrix
                backup_state(state)
                #np.savetxt(f'{out_path}/cCaps.csv', cCaps, delimiter=',')
            plt.close()
            np.savetxt(f'temp/optimized_template/matrix/matrix__{i}_{nCol}.csv', state['ks'], delimiter=',')
            #matrix = matrix + c*N
            state['ks'] = normalizeMatrix(state['ks'])
            plt.imshow(state['ks'], cmap="hot")
            plt.colorbar()
            plt.savefig(f"matrix_{i}.png")
            plt.show()
            #np.savetxt(f'{out_path}/{matrixName}_{i}_{nCol}.csv', matrix, delimiter=',')
            losses_test.append((run:=meassureBadnessOfZ(state['ks'],specs,wave,dset=photZ_test,trainOrTest="test"))[0])
            state["losses_test"] = losses_test
            plt.close()
            plt.clf()
            plt.plot(losses_train, label="train")
            plt.plot(losses_test, label="test")
            np.savetxt(f'temp/optimized_template/matrix/losses_test.csv', losses_test, delimiter=',')
            plt.legend()
            plt.savefig(f"losses.png")
            print("################### loop done #######################")
            #save plot
            #plt.savefig(f"lossByCFunc_{i}.png")
            plt.close()
            #if all ccaps are 0, break
            #state['ks'] = matrix
            if np.sum(cCaps.flatten()) == 0: 
                print("all cCaps are 0, matrix is locked!!! :O")
                break
    except KeyboardInterrupt:
        print("KeyboardInterrupt")
        pass
    print("###################DOING FINAL TEST#######################")
    losses_test.append((run:=meassureBadnessOfZ(state['ks'],specs,wave,dset=photZ_test,trainOrTest="test"))[0])
    photZ_test = run[1]
    print("################### DONE #######################")
    return state['ks'], photZ_test

