from astropy.table import Table

def load_dataframe_output(cat_out_name, template_path, out_path, out_name):
    tbl = Table.read(out_path + '/' + cat_out_name + '.zout.fits')
    restPhot = {key: tbl[key] for key in tbl.keys() if key.startswith('rest')}
    photVals = {int(key.split('rest')[-1]): tbl[key] for key in tbl.keys() if key.split('rest')[-1].isnumeric()}
    photErrs = {int(key.split('rest')[-1].split("_")[0]): tbl[key] for key in tbl.keys() if key.split('rest')[-1].split("_")[0].isnumeric() and key.endswith('err')}
    usefullColls = {
        'ID': tbl['id'],
        'z_spec': tbl['z_spec'],
        'z_phot': tbl['z_phot'],
        'z_phot_chi2': tbl['z_phot_chi2'],
        'z_phot_risk': tbl['z_phot_risk']
    }
    for filt in photVals.keys():
        usefullColls[f'{filt}_val'] = photVals[filt]
        usefullColls[f'{filt}_err'] = photErrs[filt]
    tab = Table(
        data = [d for d in usefullColls.values()],
        names = [d for d in usefullColls.keys()]
    )
    return tab

def load_dataframe_input(cat_out_name, template_path, out_path, out_name):
    tbl = Table.read(out_path + '/' + 'out_eazy' + '.fits')
    names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]
    df = tbl[names].to_pandas()
    return df

import utils_astro
import eazy_routines as ez
from IPython.utils import io

def load_photoz_output(cat_out_name, template_path, out_path, out_name, params, filts, zps, keys_id):
    # write eazy config files
    filt_num, fnames = ez.write_config(cat_out_name + '.zphot', filts, zps, keys_id,
                                        out_path=out_path, fwrite=False)

    # get a photoz object
    with io.capture_output() as captured: # capture output
        pz = ez.eazy_init_photoz(params, **fnames)
        return pz

def load_photoz_input(cat_out_name, template_path, out_path, out_name, params, filts, zps):
    # write eazy config files
    filt_num, fnames = ez.write_config(f'{cat_out_name}_{out_name}', filts, zps, [],
                                       out_path=out_path, fwrite=False)
    # get a photoz object
    with io.capture_output() as captured: # capture output
        pz = ez.eazy_init_photoz(
            params, fparam=out_path + '/' + cat_out_name + '.zphot.param',
            **fnames)
        return pz

from tqdm import tqdm

def get_output_df(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    
    df_dict = {}
    for tpath, opath, oname in tqdm(zip(templ_paths, out_paths, ftempl_strs),desc="Loading output table...", total=len(templ_paths)):
        df = load_dataframe_output(
            cat_out_name,
            template_path=tpath, out_path=opath, out_name=oname
            )
        df_dict[oname] = df
    print()
    return df_dict

def get_output_pz(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    
    pz_dict = {}
    for tpath, opath, oname in tqdm(zip(templ_paths, out_paths, ftempl_strs),desc="Loading PHOTZ output and gridding templatespace...", total=len(templ_paths)):
        params = utils_astro.gen_params(
            cat_path=cat_path, templ_path=tpath, out_path=oname,
            doCosmo=False, doUtils=False
            )
        pz = load_photoz_output(
            cat_out_name,
            template_path=tpath, out_path=opath, out_name=oname,
            params=params,
            filts=filts, zps=zps, keys_id=keys_id
            )
        pz_dict[oname] = pz
    print()
    return pz_dict

def get_input_df(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    
    df_dict = {}
    for tpath, opath, oname in tqdm(zip(templ_paths, out_paths, ftempl_strs),desc="Loading input table...", total=len(templ_paths)):
        df = load_dataframe_input(
            cat_out_name,
            template_path=tpath, out_path=opath, out_name=oname
            )
        df_dict[oname] = df
    print()
    return df_dict

def get_input_pz(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    
    pz_dict = {}
    for tpath, opath, oname, in tqdm(zip(templ_paths, out_paths, ftempl_strs),desc="Loading PHOTZ input and gridding templatespace...", total=len(templ_paths)):
        params = paramDict[oname]
        pz = load_photoz_input(
            cat_out_name,
            template_path=tpath, out_path=opath, out_name=oname,
            params=params,
            filts=filts, zps=zps
            )
        pz_dict[oname] = pz
    print()
    return pz_dict

def get_train_pz(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    
    #pz_dict = {}
    #for tpath, opath, oname, in tqdm(zip(templ_paths, out_paths, ftempl_strs),desc="Loading PHOTZ input and gridding templatespace...", total=len(templ_paths)):
    opath = "temp/eazy-output-optimizer"
    oname = "eazy"
    params = utils_astro.gen_params(
        cat_path=train_path, templ_path=matrixtemplate, out_path=opath,
        doCosmo=False, doUtils=False
        )
    filt_num, fnames = ez.write_config(f'{cat_out_name}_{oname}', filts, zps, keys_id,
                                        out_path=opath)
    pz = ez.eazy_init_photoz(
            params,
            ftran=opath + '/' + cat_out_name + "_" + oname + '.translate',
            fzp=opath + '/' + cat_out_name + "_" + oname + '.zeropoint',
            )
    #pz_dict[oname] = pz
    print()
    return pz

def get_test_pz(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, train_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate):
    return get_train_pz(templ_paths, out_paths, ftempl_strs, cat_out_name, cat_path, test_path, test_path, filts, zps, keys_id, paramDict, matrixtemplate)
    opath = "temp/eazy-output-optimizer"
    oname = "eazy"
    params = utils_astro.gen_params(
        cat_path=test_path, templ_path=matrixtemplate, out_path=oname,
        doCosmo=False, doUtils=False
        )
    filt_num, fnames = ez.write_config(f'{cat_out_name}_{oname}', filts, zps, keys_id,
                                        out_path=opath)
    pz = ez.eazy_init_photoz(
            params,
            ftran=opath + '/' + cat_out_name + "_" + oname + '.translate',
            fzp=opath + '/' + cat_out_name + "_" + oname + '.zeropoint',
            )
    print()
    return pz

