from settings import mm, fontsize, mosTiling, figwidth, DPI
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from utils_math import rebin
from matplotlib.colors import LinearSegmentedColormap

def plotSED(axis, photZs, id_cat, ftempl, ftempl_labeldict, includeChi2Val=None, logy=False,color_after_chi=False):
    global fontsize
    if len(photZs['specs'].keys()) == 0:
        print("No spectra found...")
        return
    fontsize_local = fontsize/2
    ftempl_lbl = ftempl_labeldict[ftempl]
    spec_data = photZs['specs'][id_cat]
    index = np.where(photZs['input_df'][ftempl]['id'] == id_cat)[0][0]
    pz_out = photZs['output_pz'][ftempl]
    df_out = photZs['output_df'][ftempl]
    pz_in = photZs['input_pz'][ftempl]
    df_in = photZs['input_df'][ftempl]
    chi2 = df_out['z_phot_chi2'][index]
    #chi2s = chi2s[sort]
    
    #axs = axss[i]
    #axs = axs.flatten()
    z_phot = df_out['z_phot'][df_out['ID'] == id_cat][0]
    z_spec = df_out['z_spec'][df_out['ID'] == id_cat][0]
    
    pz_out.show_fit(
        id=id_cat, add_label=False, axes=[axis],#TODO: add labels, but change config so that legend is nicer
        #show_components=True,
        zshow=z_spec,#!best guess
        xlim=[0.35, 10.0], 
        
        show_missing=True,#!these are the hexagons
        show_stars=False, snr_thresh=1.0,
        show_fnu=True,
        with_tef=True,
        #template_color='#1f77b4',
        alpha_multiplier = 1.4,
        )
    #set alphas of show_fit objects and set some markersizes
    for k,obj in enumerate(axis.lines):
        obj.set_alpha(0)
        #obj.set_linewidth(0)#!might need change
        #print color
        #print(obj.get_color())
        if k == 1: 
            obj.set_markersize(3.5)#photometry
            obj.set_alpha(0.8)
            if color_after_chi: obj.set_alpha(0)
            #remove lines arround point
            obj.set_markeredgewidth(0)
        if k == 2: 
            obj.set_markersize(3.5)#photometry, bad points?
            obj.set_alpha(0.25)
            obj.set_markeredgewidth(0)
        if k == 0: 
            #obj.set_marker('x')
            obj.set_alpha(0)
        if k == 3: 
            obj.set_linewidth(1)#best fit
            obj.set_alpha(0.5)
            obj.set_zorder(19)
            #if color_after_chi: obj.set_alpha(0)
        #obj.set_visible(False)
    #and errorbars
    for k,obj in enumerate(axis.collections):
        obj.set_alpha(0.5)#errorbar
        #obj.set_linewidth(0)#!might need change
        if color_after_chi: obj.set_alpha(0)
        """if k == 5:
            obj._edgecolors = np.zeros(obj._edgecolors.shape)
            #set color pink
            obj.set_facecolor('#ff7f0e')
        obj._edgecolors = np.zeros(obj._edgecolors.shape)
        obj.set_facecolor('#ff7f0e')"""
        if k not in [1]:
            obj.set_visible(False)
            pass
        else:
            pass
            size = obj.get_sizes()
            obj.set_sizes(size*0.5)
            obj.set_linewidth(0)
        if k in [3]:
            obj.set_visible(True)
    
    data_phot = pz_out.show_fit(id=id_cat, zshow=z_phot, get_spec=True, show_fnu=True)
    #fit_chi2red = data_phot['chi2']/(sum(data_phot['valid'])-1)
    

        
    #and fill_between (polygons)
    for obj in axis.patches:
        obj.set_alpha(0)#!temp
        obj.set_linewidth(0)#!might need change

    
    #!plot actuall guess
    #TODO: fix
    
    tempspec_zphot = [np.array((data_phot['templz']*u.AA).to(u.um)), np.array((data_phot['templf']*u.uJy).to(u.uJy))]
    #limit to axis.get_xlim()
    tempspec_zphot = [tempspec_zphot[0][tempspec_zphot[0] > axis.get_xlim()[0]], tempspec_zphot[1][tempspec_zphot[0] > axis.get_xlim()[0]]]
    tempspec_zphot = [tempspec_zphot[0][tempspec_zphot[0] < axis.get_xlim()[1]], tempspec_zphot[1][tempspec_zphot[0] < axis.get_xlim()[1]]]
    #rebin to simular wavelength intervals
    intervals = np.array([tempspec_zphot[0][k+1] - tempspec_zphot[0][k] for k in range(len(tempspec_zphot[0])-1)])
    meanInterval = np.mean(intervals)
    targetWave = copy(spec_data['wave'])
    targetWave = np.linspace(targetWave[0], targetWave[-1], len(targetWave))
    #use rebin function
    #tempspec_zphot = [targetWave, np.interp(targetWave, tempspec_zphot[0], tempspec_zphot[1])]
    flux_rebin, _ = rebin(tempspec_zphot[0], tempspec_zphot[1], targetWave, np.ones(len(data_phot['templf']), dtype=np.float64))#!Don't think there is a error for best template fit? Please check:)
    tempspec_zphot = [targetWave, flux_rebin]
    axis.plot(tempspec_zphot[0], tempspec_zphot[1], c='tab:orange', lw=1, ls='--',alpha=0.9,zorder=20)
    #if i == 0: axis.set_title(f'{includeChi2s[j]*100:.0f}%-tile $\chi^2$', fontsize=fontsize*2)
    if includeChi2Val != None:
        plt.title(f'{includeChi2Val*100:.0f}%-tile $\chi^2$', fontsize=fontsize_local*2)
    #set y max lim to 1.5 times higher
    

    
    
    #annotate to the right of the frame; chi2, objid, redshift, delta_z
    anXs = [1.02]*5
    anYs = np.linspace(0.0, len(anXs)*0.07, len(anXs))+(0.5-(len(anXs)*0.07/2))
    chi2red = chi2/(sum(data_phot['valid'])-1)
    modelChi2red = np.sum((data_phot['model'][data_phot['valid']]-data_phot['fobs'][data_phot['valid']])**2/(data_phot['efobs'][data_phot['valid']]**2))/(sum(data_phot['valid'])-1)
    print("Chi2 diff:", chi2red-modelChi2red)
    axis.annotate(f'{ftempl_lbl}', xy=(anXs[4],anYs[4]), xycoords='axes fraction', fontsize=fontsize_local*1.25, ha='left', va='center', textcoords='offset points', xytext=(0,0))#TODO: move to left
    axis.annotate(f'$\chi^2_{{red}}$≈{"{:.1e}".format(chi2red)}', xy=(anXs[0],anYs[0]), xycoords='axes fraction', fontsize=fontsize_local, ha='left', va='center', textcoords='offset points', xytext=(0,0))
    axis.annotate(f'ID={id_cat}', xy=(anXs[1],anYs[1]), xycoords='axes fraction', fontsize=fontsize_local*0.8, ha='left', va='center', textcoords='offset points', xytext=(0,0))
    axis.annotate(f'$z_{{spec}}$≈{z_spec:.1f}', xy=(anXs[2],anYs[2]), xycoords='axes fraction', fontsize=fontsize_local, ha='left', va='center', textcoords='offset points', xytext=(0,0))
    axis.annotate(f'$\Delta z$≈{z_phot-z_spec:.1f}', xy=(anXs[3],anYs[3]), xycoords='axes fraction', fontsize=fontsize_local, ha='left', va='center', textcoords='offset points', xytext=(0,0))

    #plot spec data
    
    axis.plot(spec_data['wave'], spec_data['flux'], c='k', lw=0.5, ls='-', alpha=0.6, label="Current spectra")#!added this, validate
    axis.fill_between(spec_data['wave'], spec_data['flux']-spec_data['flux_err'], spec_data['flux']+spec_data['flux_err'], color='k', alpha=0.2, lw=0, zorder=25)

    #rebin model spectra to simular wavelength intervals
    intervals = np.array([spec_data['wave'][k+1] - spec_data['wave'][k] for k in range(len(spec_data['wave'])-1)])
    meanInterval = np.mean(intervals)
    

    x_fit = axis.lines[3].get_xdata()#where best fit template is stored
    y_fit = axis.lines[3].get_ydata()
    try:
        x_fit = x_fit.to(u.AA).value
        y_fit = y_fit.to(u.uJy).value
        print("Found odd unit thing:P")
    except:
        pass
    x_fill = axis.collections[5].get_paths()[0].vertices[:,0]#where best fit template error is stored
    y_fill = axis.collections[5].get_paths()[0].vertices[:,1]
    indMaxxfill = np.where(x_fill == x_fill.max())[0][0]
    xfillsort = np.argsort(x_fill)
    x_fill = x_fill[xfillsort]
    y_fill = y_fill[xfillsort]
    x_fill_bot = []
    y_fill_bot = []
    x_fill_top = []
    y_fill_top = []
    for k in range(len(x_fill)):#lookup by uniqueness
        ys = []
        x = x_fill[k]
        if x not in x_fill_bot:
            x_fill_bot.append(x)
            x_fill_top.append(x)
            for l in range(k, len(x_fill)):
                if x_fill[l] == x:
                    ys.append(y_fill[l])
                else:
                    break
            #y_fill_bot.append(np.min(ys))
            """if len(ys)//2 != len(ys)/2:
                raise ValueError("Unexpected format of y_fill")"""
            """for j in range(0, len(ys),2):
                y_fill_top.append(max(ys[j], ys[j+1]))
                y_fill_bot.append(min(ys[j], ys[j+1]))
                if j > 0:
                    x_fill_top.append(x)
                    x_fill_bot.append(x)"""
            if len(ys) <= 1: raise ValueError("Unexpected format of y_fill")
            y_fill_top.append(np.max(ys))
            y_fill_bot.append(np.min(ys))
            
    x_fill_bot = np.array(x_fill_bot)
    y_fill_bot = np.array(y_fill_bot)
    x_fill_top = np.array(x_fill_top)
    y_fill_top = np.array(y_fill_top)

    y_fill_delta = (y_fill_top - y_fill_bot)/2
    xMax = axis.get_xlim()[1]
    xMin = axis.get_xlim()[0]
    axis.set_xlim((xMin:=xMin), (xMax:=5))#!chaning xMin ruins the xticks
    axis.set_ylim(-axis.get_ylim()[1]*0.1, axis.get_ylim()[1]*1)
    topSort = np.argsort(x_fill_top)
    botSort = np.argsort(x_fill_bot)
    x_fill_top = x_fill_top[topSort]
    y_fill_top = y_fill_top[topSort]
    x_fill_bot = x_fill_bot[botSort]
    y_fill_bot = y_fill_bot[botSort]
    """plt.plot(x_fill_top, y_fill_top)
    plt.plot(x_fill_bot, y_fill_bot)
    plt.xlim(xMin, xMax)
    plt.ylim(*axis.get_ylim())
    plt.show()"""

    #rebin to target intervals
    targetFlux, _ = rebin(x_fit, y_fit, targetWave,y_fill_delta)
    targetFlux_yfill_top = np.interp(targetWave, x_fill_top, y_fill_top)
    targetFlux_yfill_bot = np.interp(targetWave, x_fill_bot, y_fill_bot)
    """plt.plot(targetWave, targetFlux_yfill_top)
    plt.plot(targetWave, targetFlux_yfill_bot)
    plt.show()"""

    addWave = targetWave[-1] + meanInterval
    while addWave < x_fit[-1] and addWave < xMax:
        targetWave = np.append(targetWave, addWave)
        targetFlux = np.append(targetFlux, np.interp(addWave, x_fit, y_fit))
        targetFlux_yfill_top = np.append(targetFlux_yfill_top, np.interp(addWave, x_fill_top, y_fill_top))
        targetFlux_yfill_bot = np.append(targetFlux_yfill_bot, np.interp(addWave, x_fill_bot, y_fill_bot))
        addWave += meanInterval
    addWave = targetWave[0] - meanInterval
    while addWave > x_fit[0] and addWave > xMin:
        targetWave = np.insert(targetWave, 0, addWave)
        targetFlux = np.insert(targetFlux, 0, np.interp(addWave, x_fit, y_fit))
        targetFlux_yfill_top = np.insert(targetFlux_yfill_top, 0, np.interp(addWave, x_fill_top, y_fill_top))
        targetFlux_yfill_bot = np.insert(targetFlux_yfill_bot, 0, np.interp(addWave, x_fill_bot, y_fill_bot))
        addWave -= meanInterval
    """plt.plot(targetWave, targetFlux_yfill_top)
    plt.plot(targetWave, targetFlux_yfill_bot)
    plt.show()"""
    #stitch errorverticies together again
    targetFlux_xfill = np.append(targetWave, targetWave[::-1])
    targetFlux_yfill = np.append(targetFlux_yfill_bot, targetFlux_yfill_top[::-1])

    axis.lines[3].set_xdata(targetWave)
    axis.lines[3].set_ydata(targetFlux)
    axis.collections[5]._paths[0].verticies = np.array([targetFlux_xfill, targetFlux_yfill]).T

    #if y is badly scaled:
    ylim_max = axis.get_ylim()[1]
    ylim_min = axis.get_ylim()[0]
    xlim_max = axis.get_xlim()[1]
    #see if ylim_max is higher than targetFlux within panel
    if ylim_max > np.max(targetFlux[np.where(targetWave <  xlim_max)]):
        axis.set_ylim(ylim_min, np.max(targetFlux[np.where(targetWave <  xlim_max)]))

    #if ylog
    if logy:
        axis.set_yscale('log')
        axis.set_ylim(ylim_max/100, ylim_max*5)

    data_spec = pz_out.show_fit(id=id_cat, zshow=z_spec, get_spec=True, show_fnu=True)
    if color_after_chi:
        #replace the photometry points with a color representing that points chi
        #get chis
        v = data_spec['valid']
        fobs = data_spec['fobs'][v]
        ferr = data_spec['efobs'][v]
        fgues = data_spec['model'][v]
        w = data_spec['pivot'][v]
        mask = np.where(fobs >= axis.get_ylim()[0]) or np.where(fgues >= axis.get_ylim()[0])
        fobs = fobs[mask]
        ferr = ferr[mask]
        fgues = fgues[mask]
        w = w[mask]
        chi2s = [(p-g)**2/(e**2) for p,g,e in zip(fobs,fgues,ferr)]
        chi2s = np.array(chi2s)
        chi2s = chi2s[~np.isnan(chi2s)]
        print(np.sum(chi2s))
        chi2s = np.log10(chi2s)+1
        print(np.sum(chi2s)/len(chi2s))
        colors = [(0, 0, 0), (1, 0, 0)] # first color is black, last is red
        cmap = LinearSegmentedColormap.from_list(
            "Custom", colors, N=20)
        for i, chi2 in enumerate(chi2s):
            axis.scatter(w[i]/10000, fobs[i], c=cmap(chi2/chi2s.max()), s=10, zorder=20)
            axis.scatter(w[i]/10000, fgues[i], c='tab:orange', s=10, zorder=19)
            axis.errorbar(w[i]/10000, fobs[i], yerr=ferr[i], c=cmap(chi2/chi2s.max()), zorder=20)

    #set non-log axis
    axis.set_xscale('linear')
    ylabel = axis.get_ylabel().replace('\\mu', 'n')
    axis.set_ylabel(ylabel, fontsize=fontsize_local)
    axis.set_yticklabels([f'{tick*1e3:.0f}' for tick in axis.get_yticks()], fontsize=fontsize_local)
    axis.set_xticklabels([f'{tick:.1f}' for tick in axis.get_xticks()], fontsize=fontsize_local)
    axis.set_xlabel('$\lambda_{{obs}}$ [$\mu$m]', fontsize=fontsize_local)
    
from math import ceil
from copy import copy

def plot_SED_mosaic(photZs,ftempl_labels,ftempl_strs,ftempl_labeldict,runTime=0,idx=None):
    #fig, axss = plt.subplots(len(ftempl_strs), 2, figsize=(16, 4*len(ftempl_strs)), dpi=300, width_ratios=[2, 1])
    if idx == None:
        includeChi2s = [0.1,0.5,0.8,0.9,1]
    else: includeChi2s = np.full(len(idx), np.nan)
    
    figwidth = 183*mm
    global fontsize
    fontsize_local = fontsize/2
    plotSucces = False

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

        figs, axss = [], []
        for chi in includeChi2s:
            figMos, axsMos = plt.subplots(*tiling, figsize=figsize, dpi=DPI, facecolor=(1,1,1,0))
            if mosaicLen == 1: axsMos = [axsMos]
            if len(mode_ftempl_lbls) == 1: axsMos = [axsMos]
            figs.append(figMos)
            axss.append(axsMos)
            plt.subplots_adjust(hspace=0, wspace=0.15*(len(includeChi2s)-1))
            plt.close(figMos)
    
        mode_output_df = [photZs['output_df'][ftempl] for ftempl in mode_ftempl_strs]    
        for i, df_out, ftempl in zip(range(len(mode_ftempl_strs)), mode_output_df, mode_ftempl_strs):
            #print(ftempl)
            if len(photZs['specs'].keys()) == 0: 
                print("No spectra found...")
                continue
            if idx == None:#have sort by chi2
                chi2s = df_out['z_phot_chi2']
                chi2s = [chi2 if chi2 > 0 else 0 for chi2 in chi2s]
                chi2s = [chi2 if df_out['z_spec'][i] > 0 else 0 for i, chi2 in enumerate(chi2s)]
                ids = df_out['ID']
                chi2s = [chi2 if id in photZs['specs'].keys() else 0 for id, chi2 in zip(ids, chi2s)]
                sort = list(np.argsort(chi2s))[chi2s.count(0):]
                #remove if chi2 < 0:
                include = [int(inc*(len(sort)-1)) for inc in includeChi2s]
                sort = np.array(sort)[include]
            else:
                include = []
                for inc in idx:
                    if len((index:=np.where(df_out['ID'] == inc)[0])) > 0:
                        include.append(index[0])
                sort = np.array(include)
            del include
            for j in range(sort.shape[0]-1,0-1,-1):
                #if spec not exist
                if df_out['ID'][sort[j]] not in photZs['specs'].keys():
                    sort = np.delete(sort, j)
            if len(sort) == 0: 
                print(f"no spec for {ftempl}")
                continue
            IDs = df_out['ID'][sort]
            if idx == None: chi2s = np.array(chi2s)[sort]
            else: chi2s = np.full(len(sort), np.nan)
            

            for j, chi2, id_cat in zip(range(len(sort)), chi2s, IDs):
                #print(f"Chi2: {chi2:.2f} ID: {id_cat}")
                axis = axss[j][i//int(mosTiling)][i%int(mosTiling)]
                
                plotSED(axis, photZs, id_cat, ftempl, ftempl_labeldict, includeChi2s[j])
                if i//mosTiling != len(axss[j])-1: 
                    axis.set_xlabel('')
                    axis.set_xticklabels([])
                else: axis.set_xlabel('$\lambda$ [$\mu$m]', fontsize=fontsize_local)
            plotSucces = True
                
        
        if not plotSucces: continue
        plt.clf()
        plt.show()
        for i,fig,chi2 in zip(range(len(figs)),figs,includeChi2s):
            if idx != None: chi2 = idx[i]
            fig.savefig(f'./figures/forpaper/seds_{mode}_{chi2}_{runTime}.png', dpi=300, bbox_inches='tight', transparent=True)
            fig.clf()
            plt.close(fig)
            #show by loading image
            img = plt.imread(f'./figures/forpaper/seds_{mode}_{chi2}_{runTime}.png')
            fig = plt.figure(figsize=(16, 4*len(ftempl_strs)), dpi=100, facecolor=(1,1,1,0))
            plt.axis('off')
            plt.imshow(img)
            plt.show()