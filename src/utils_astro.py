

####################### DATA #######################

filts_nircam = {
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
}

ftempl_labeldict = {#! use ; to indicate modifications
    'blue_sfhz_13': "BLSFH",
    'carnall_sfhz_13': "CASFH",
    'corr_sfhz_13': "COSFH",
    'eazy_v1.3.spectra': "EAZ3",
    'EMextreme': "EMEx",
    'EMlines': "EMLi",
    'EMLines;linearcomb': "EMLi;LC",
    'fsps_45k': "F45k",
    'fsps_60k': "F60k",
    'fsps_45k;linearcomb': "F45k;LC",
    'fsps_45k;0.3removed': "F45k;0.3r",
}

"""catalogue_paths = {
    "hlsp_jades_jwst_nircam_goods-s-deep_photometry_v1.0_catalog_large_withSpec.fits": "gds/jades/initial/phot"

}"""

####################### FUNCTIONS #######################

def get_nircam_filters():
    return filts_nircam

def get_HST_filters():
    return filts_HST

def get_all_filters():
    filts_HST = get_HST_filters()
    filts_nircam = get_nircam_filters()
    filts = {**filts_nircam, **filts_HST}
    return filts

from eazy import filters
def get_filts_wavelengths_from_FILTERRES():
    filt_wave = {}
    res = filters.FilterFile('FILTER.RES.latest')
    for f in res.filters:
        for filt in filts_nircam | filts_HST:
            if filt.lower() in f.name:
                filt_wave[filt] = {
                    'mean': f.pivot,
                    'width': f.rectwidth
                    }
                break
    return filt_wave

import os
import numpy

def get_templates():
    ftempl_strs = []
    ftempl_labels = []
    files = numpy.sort(os.listdir("templates2test"))
    for f in files:
        if f.endswith(".param"):
            templ_str = f.split(".param")[0]
            if templ_str not in ftempl_labeldict.keys(): continue
            ftempl_strs.append(templ_str)
            try:
                ftempl_labels.append(ftempl_labeldict[templ_str])
            except KeyError:
                raise KeyError("No label for template: ", templ_str)
        else: continue
        print("Found template: ", templ_str)
        print("Label: ", ftempl_labeldict[templ_str])
        print("")
    templ_paths = [f"templates/{e}.param" for e in ftempl_strs]
    return ftempl_strs, ftempl_labels, ftempl_labeldict, templ_paths

import helper_module as hmod
import eazy_routines as ez
import numpy as np
from astropy import units as u
from astropy.table import join as jointab
from astropy.table import Table

def get_mv_reddening():
    filts = get_all_filters()
    mw_reddening = ez.get_atten_dict(filts)
    return mw_reddening

def get_hdu_names(catalogue_inpath):
    hdu_names = []
    with fits.open(catalogue_inpath) as hdul:
        for h in hdul:
            try:
                hdu_names.append(h.name)
            except AttributeError:
                hdu_names.append(None)
    return hdu_names

def get_secondarry_zpec(reverting_catalogues):
    print("ATTENTION!: No z_spec column found in catalogue. Finding z_spec from catalogue alternatives.")
    for cat in reverting_catalogues:
        hdu_names = get_hdu_names(cat)
        for i,h in enumerate(hdu_names):
            if i == 0: continue
            print("Trying hdu: ", h, ", in catalogue: ", cat)
            tab = Table.read(cat, hdu=i)
            cols = [str(c).lower() for c in tab.colnames]
            if 'z_spec' in cols:
                print("Found z_spec in: ", cat)
                zspec_naming = tab.colnames[np.where(np.array(cols) == 'z_spec')[0][0]]
                id_namings = ["NIRCam_ID", "nircam_id", "id", "ID"]
                for id_naming in id_namings:
                    if id_naming in tab.colnames:
                        #rename naming to ID, z_spec
                        tab.rename_column(id_naming, 'ID')
                        tab.rename_column(zspec_naming, 'z_spec')
                        #remove rows with masked "ID" or "z_spec"
                        try: 
                            tab = tab[~tab['ID'].mask]
                            tab = tab[~tab['z_spec'].mask]
                        except AttributeError: pass
                        return tab['ID', 'z_spec']
    print("No z_spec found in any catalogues!")
    raise ValueError("No z_spec found in any catalogues!")
    
def catalogue_2_eazytable(catalogue_inpath:str, cat_out_name, reverting_catalogues=[], z_min_limit=0.0):
    mw_reddening = get_mv_reddening()
    hdu_names = get_hdu_names(catalogue_inpath)
    hdu_tab_in = np.where(np.array(hdu_names) == 'CIRC')[0][0]
    hdu_tab_redshifts = np.where(np.array(hdu_names) == 'PHOTOZ')[0][0]
    tab_in = Table.read(catalogue_inpath, hdu=hdu_tab_in)
    tab_redshifts = Table.read(catalogue_inpath, hdu=hdu_tab_redshifts)
    # load fluxes
    # CIRC1: 0.10 arcsec aperture (see README)
    ext = '_CIRC1'
    cols_dummy = hmod.get_matches(ext, tab_in.columns, exclude='_ei')
    cols_f = np.sort(hmod.get_matches(ext, cols_dummy, exclude='_e'))
    cols_fe = np.sort(hmod.get_matches('_e', cols_dummy))
    cols_fluxes = list(np.vstack([cols_f, cols_fe]).T.flatten())
    cols = list(np.insert(cols_fluxes, 0, ['ID', 'RA', 'DEC', 'z_spec']))

    #check wether z_spec is available
    if 'z_spec' not in tab_redshifts.colnames:
        tab_redshifts = get_secondarry_zpec(reverting_catalogues)
    #print len of z_spec
    print("z_spec length: ", len(tab_redshifts['z_spec']))
    tab_in = jointab(tab_in, tab_redshifts['ID', 'z_spec'], join_type='inner', keys='ID')
    tab_out = tab_in[cols]#! Is this filtering fluxes out of filters that might have data??

    # convert from nJy to uJy
    # and apply MW reddening
    keys = np.array(list(mw_reddening.keys()))
    for c in cols_fluxes:
        # convert from nJy to uJy
        tab_out[c].unit = u.nJy
        tab_out[c] = tab_out[c].to(u.uJy)
        
        # apply MW reddening
        matches = hmod.get_matches(keys, c, get_idxs=True)
        key = keys[np.int32(matches[:,0])][0]
        tab_out[c] *= mw_reddening[key]

    # redshift limit

    for i in range(len(tab_out)):
        if tab_out['z_spec'][i] < z_min_limit:
            tab_out['z_spec'][i] = -1.0
        
    # rename columns
    for c in cols_f:
        cnew = c.replace(ext, '_flux')
        tab_out.rename_column(c, cnew)

    for c in cols_fe:
        cnew = c.replace(ext+'_e', '_err')
        tab_out.rename_column(c, cnew)

    #=== apply MW reddening
    #atten_dict = ez.get_atten_dict(filts_eazyres, filts_str)
    #degr_image_sig *= atten_dict[filt] / 100. # uJy

    # save EAZY table
    tab_out.write(f'temp/{cat_out_name}.fits', format='fits', overwrite=True)

    #make a keys_id
    keys_id = ['ID id', 'RA ra', 'DEC dec', 'z_spec z_spec']

    return f"temp/{cat_out_name}.fits", keys_id

from astropy.cosmology import Planck18
def gen_params(cat_path,templ_path,out_path,cosmo=Planck18,maxZ=12, doCosmo=True, doUtils=True):
    params = {
        "cat_path": cat_path,
        "templ_path": templ_path,
        "out_path": out_path,
        "FIX_ZSPEC": 'n'
        }
    if doCosmo:
        addon = {
            "Z_MAX": 15.0,
            "H0": cosmo.H0,
            "OMEGA_M": cosmo.Om0,
            "OMEGA_L": cosmo.Ode0
        }
    else:
        addon = {}
    params = {**params, **addon}
    if doUtils:
        addon = {
            "CATALOG_FORMAT": 'fits',
            'VERBOSE': 'y',
            "USE_ZSPEC_FOR_REST": 'n'
        }
    else:
        addon = {}
    params = {**params, **addon}
    return params

def get_outpaths(eazy_out, cat_out_name, ftempl_strs, runtimeNum=-1):
    runTimes = np.sort([int(f.split('_')[-1]) for f in os.listdir(eazy_out) if f != '.gitignore'])
    runTimes = np.unique(runTimes)
    runTime = int(runTimes[runtimeNum])
    print("Picking runTime:", runTime)
    outpaths_f = '{eazy_outfolder}/{ftempl}_{runTime}/' + cat_out_name + '.zout.fits'
    out_paths = [os.path.dirname(outpaths_f.format(ftempl=f, runTime=runTime, eazy_outfolder=eazy_out)) for f in ftempl_strs]
    return out_paths

import pandas as pd
from astropy.io import fits
from copy import copy

def get_spectra(spec_dir, input_df, ftempl_strs):
    fpath = os.path.join(os.getenv('astrodata'), spec_dir)
    print(os.listdir(fpath))
    fnames_1d = [f for f in os.listdir(fpath) if f.endswith('1d.fits')]
    fpaths_1d = np.sort([os.path.join(fpath, f) for f in fnames_1d])

    df_spec = pd.DataFrame(columns=['ID', 'RA', 'DEC'], dtype=np.float32)
    spec_data = {}
    for f in fpaths_1d:
        id_fname = int(f.split('/')[-1].split('-')[3].strip('_clear'))
        with fits.open(f) as hdul:
            hdr = hdul[0].header
            row = np.array([id_fname, hdr['RA'], hdr['DEC']]).reshape(1, -1)
            df_cur = pd.DataFrame(row, columns=['ID', 'RA', 'DEC'])
            df_spec = pd.concat([df_spec, df_cur])
            
            
            unit = u.Unit(hdr['BUNIT']) # erg/s/cm2/A
            flux = hdul[1].data['FLUX'] * unit
            fluxerr = hdul[1].data['FLUX_ERR'] * unit
            wave = hdul[1].data['WAVELENGTH'] * u.um
            equiv = u.spectral_density(wave.to(u.AA))
            flux = flux.to(unit * u.AA / u.Hz, equivalencies=equiv)
            fluxerr = fluxerr.to(unit * u.AA / u.Hz, equivalencies=equiv)
            spec_data[id_fname] = {
                'wave': wave.value, 
                'flux': flux.to(u.uJy).value, 
                'flux_err': fluxerr.to(u.uJy).value
            }
            #spec_data[id_fname] = hdul[1].data['FLUX']
        
    df_spec.reset_index(drop=True, inplace=True)
    df_spec.ID = df_spec.ID.astype(np.int32)

    # match spec sample to photo catalog
    from astropy.coordinates import SkyCoord
    def match_catalogs(samp_x, samp_y, cat_x, cat_y, max_sep=1.0*u.arcsec):
        sample = SkyCoord(ra=samp_x*u.degree, dec=samp_y*u.degree)
        catalog = SkyCoord(ra=cat_x*u.degree, dec=cat_y*u.degree)
        idx, d2d, d3d = sample.match_to_catalog_sky(catalog)
        samp_sel = d2d < max_sep
        return samp_sel, idx

    #photZs['input_df'][ftempl_strs[0]] = photZs['input_df'][ftempl_strs[0]]
    cat_x = input_df[ftempl_strs[0]]['ra'].values
    cat_y = input_df[ftempl_strs[0]]['dec'].values
    samp_x = df_spec['RA'].values
    samp_y = df_spec['DEC'].values
    mask_samp, idx_cat = match_catalogs(samp_x, samp_y, cat_x, cat_y, 
                                        max_sep=0.8*u.arcsec)

    # add spec-flag column to catalog
    input_df_updated = copy(input_df)
    for i in range(len(ftempl_strs)):
        input_df_updated[ftempl_strs[i]]['ID_spec'] = np.full(len(input_df_updated[ftempl_strs[0]]), -1, dtype=np.int32)
        input_df_updated[ftempl_strs[i]].loc[idx_cat[mask_samp], 'ID_spec'] = df_spec.loc[mask_samp, 'ID'].values

    #refrhase spec_data keys as phot ID instead
    temp = {}
    for i,k in enumerate(spec_data):
        if mask_samp[i] == False: continue
        spec_data[k]['ID_spec'] = k
        indexPhot = idx_cat[i]
        photID = input_df_updated[ftempl_strs[0]]['id'].values[indexPhot]
        temp[photID] = spec_data[k]
    spec_data = temp

    return spec_data, input_df_updated

import utils_math

def clean_spectra(wave, flux, flux_err):
    #find nans in flux and flux_err
    nan_index = np.where(np.isnan(flux) | np.isnan(flux_err) | np.isnan(wave))
    #remove nans
    wave = np.delete(wave, nan_index)
    flux = np.delete(flux, nan_index)
    flux_err = np.delete(flux_err, nan_index)
    return wave, flux, flux_err

def rebin_spectra(spec_data, rebinWidth=0.1):
    
    wave = spec_data['wave']
    flux = spec_data['flux']
    flux_err = spec_data['flux_err']
    #remove nans
    """for i in range(len(wave)-1,-1,-1):
        if np.isnan(flux[i]) or np.isnan(flux_err[i]):
            wave = np.delete(wave, i)
            flux = np.delete(flux, i)
            flux_err = np.delete(flux_err, i)"""
    wave, flux, flux_err = clean_spectra(wave, flux, flux_err)
    #wave_new, flux_new, flux_err_new = rebin_spec(wave, flux, flux_err, xlo=wave.min(), xhi=wave.max(), bw=rebinWidth)
    wave_new = np.arange(wave.min(), wave.max()+rebinWidth, rebinWidth)
    flux_new, flux_err_new = utils_math.rebin(wave, flux, wave_new, y_err=flux_err)
    spec_data['wave'] = wave_new
    spec_data['flux'] = flux_new
    spec_data['flux_err'] = flux_err_new  
    return spec_data

def scale_spec_to_phot(spec_data,phot_df,phot_pz,IDkey,ftempl_strs,filts):
    wave = spec_data['wave']*u.um
    flux = spec_data['flux']
    flux_err = spec_data['flux_err']
    #get photometry by fnu in photoz object
    z_spec = phot_df[ftempl_strs[0]]['z_spec'].values
    index = np.where(phot_df[ftempl_strs[0]]['id'] == IDkey)[0][0]
    data = phot_pz[ftempl_strs[0]].show_fit(id=IDkey, zshow=z_spec[index], show_fnu=True, get_spec=True)
    fobs = data['fobs']*u.uJy
    
    filts_dict = phot_pz[ftempl_strs[0]].filters
    selection = list(range(len(filts)))#[7,8,9,10,11,12,13,14]#only selects the middle bunch of filters
    filts_dict = [filts_dict[i] for i in range(len(filts_dict)) if i in selection]
    fobs = [fobs[i] for i in range(len(fobs)) if i in selection]
    wave = np.array([wave[i].to(u.AA).value for i in range(len(wave))])
    fobs = [fobs[i].to(u.uJy).value for i in range(len(fobs))]
    for i in range(len(fobs)-1,-1,-1):
        if fobs[i] < -10:
            fobs = np.delete(fobs, i)
            filts_dict = np.delete(filts_dict, i)

    #find scalar scaling
    filtered_spectravals = []
    for i,filt in enumerate(filts_dict):
        wave_filt = (filt.wave*u.AA).to(u.AA).value
        wave_throughput = filt.throughput
        filterValue = utils_math.bandpass(wave,flux,wave_filt,wave_throughput)
        filtered_spectravals.append(filterValue)
        #print(filterValue)
    #normalize spectra
    fnus_procc = copy(fobs)
    for i in range(len(fobs)-1,-1,-1):
        if np.isnan(filtered_spectravals[i]):
            fnus_procc = np.delete(fnus_procc, i)
            filtered_spectravals = np.delete(filtered_spectravals, i)
    if np.sum(filtered_spectravals) == 0: raise ValueError("SUS!")
    scaling = (np.sum(fnus_procc) / np.sum(filtered_spectravals))
    flux *= scaling
    flux_err *= scaling
    
    #save
    spec_data['flux'] = flux
    spec_data['flux_err'] = flux_err
    return spec_data