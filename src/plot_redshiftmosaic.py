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
zmin = zmin
zMax = zMax
mosTiling = mosTiling
mm = eval(mm)
figwidth = eval(figwidth)
DPI = DPI
ro = ro
zCharacteristic = zCharacteristic







#TODO: needs tiling color by accumulative flux

etasStartIndex = 0
etas = []

from copy import copy, deepcopy
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

def plot(df_out,ftempl_strs,ftempl_labels,runTime):
    global etasStartIndex
    global etas
    etasStartIndex = 0
    etas = []
    #precalc stats
    deltaZ_s = {'default':[], 'modified':[]}
    deltaZ_all_s = {'default':[], 'modified':[]}
    outliers_s = {'default':[], 'modified':[]}
    outliers_all_s = {'default':[], 'modified':[]}
    bias_s = {'default':[], 'modified':[]}
    scatter_s = {'default':[], 'modified':[]}
    etas_s = {'default':[], 'modified':[]}
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
            redshiftTbl = df_out[ftempl]
            redChar = {
                'z_spec': np.array(redshiftTbl['z_spec']),
                'z_phot': np.array(redshiftTbl['z_phot']),
            }
            mask = (redChar['z_spec'] > zCharacteristic) & (~np.isnan(redChar['z_spec'])) & (~np.isnan(redChar['z_phot']) & (redChar['z_phot'] > 0))
            redChar['z_spec'] = redChar['z_spec'][mask]
            redChar['z_phot'] = redChar['z_phot'][mask]
            deltaZ = (redChar['z_phot'] - redChar['z_spec'])/(1 + redChar['z_spec'])
            deltaZ_s[mode].append(deltaZ)
            deltaZ_all = (redshiftTbl['z_phot'] - redshiftTbl['z_spec'])/(1 + redshiftTbl['z_spec'])
            deltaZ_all_s[mode].append(deltaZ_all)
            outliers = (deltaZ > 0.15) | (deltaZ < -0.15)
            outliers_s[mode].append(outliers)
            outliers_all = (deltaZ_all > 0.15) | (deltaZ_all < -0.15)
            outliers_all_s[mode].append(outliers_all)
            bias = np.mean(deltaZ)
            bias_s[mode].append(bias)
            scatter = np.std(deltaZ)
            scatter_s[mode].append(scatter)
            eta = np.sum(outliers)/len(outliers)
            etas_s[mode].append(eta)

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
        else:
            axsMos = [axsMos]
        plt.subplots_adjust(wspace=0, hspace=0)
        
        """#precalc stats
        deltaZ_s = []
        deltaZ_all_s = []
        outliers_s = []
        outliers_all_s = []
        bias_s = []
        scatter_s = []
        etas_s = []
        for i, ftempl in enumerate(mode_ftempl_strs):
            if ftempl not in photZs['output_df'].keys():
                continue
            redshiftTbl = photZs['output_df'][ftempl]
            redChar = {
                'z_spec': np.array(redshiftTbl['z_spec']),
                'z_phot': np.array(redshiftTbl['z_phot']),
            }
            mask = (redChar['z_spec'] > zCharacteristic) & (~np.isnan(redChar['z_spec'])) & (~np.isnan(redChar['z_phot']) & (redChar['z_phot'] > 0))
            redChar['z_spec'] = redChar['z_spec'][mask]
            redChar['z_phot'] = redChar['z_phot'][mask]
            deltaZ_s.append((redChar['z_phot'] - redChar['z_spec'])/(1 + redChar['z_spec']))
            deltaZ = deltaZ_s[-1]
            deltaZ_all_s.append((redshiftTbl['z_phot'] - redshiftTbl['z_spec'])/(1 + redshiftTbl['z_spec']))
            deltaZ_all = deltaZ_all_s[-1]
            outliers_s.append((deltaZ > 0.15) | (deltaZ < -0.15))
            outliers = outliers_s[-1]
            outliers_all_s.append((deltaZ_all > 0.15) | (deltaZ_all < -0.15))
            outliers_all = outliers_all_s[-1]
            bias_s.append(np.mean(deltaZ))
            bias = bias_s[-1]
            scatter_s.append(np.std(deltaZ))
            scatter = scatter_s[-1]
            etas_s.append(np.sum(outliers_s[i])/len(outliers_s[i]))
            eta = etas_s[-1]"""

        #color the backgrounds of the plots with colormap RdYlGn by eta
        temp = [etas_s[k] for k in etas_s.keys()]
        etas = []
        for i in range(len(temp)): etas += temp[i]
        
        minEr, maxEr = np.min(etas)**(-ro), np.max(etas)**(-ro)
        norm = mpl.colors.Normalize(vmin=maxEr, vmax=minEr)
        cmap = mpl.cm.RdYlGn
        cmap.set_under(color=cmap(maxEr))
        for i in range(etasStartIndex,len(axsMos)+etasStartIndex):
            j = i - etasStartIndex
            color = cmap(norm(etas[i]**(-ro)))
            color = (*color[:-1], 0.1)#set alpha = 0.1
            axsMos[j].set_facecolor(color)
        etasStartIndex = len(axsMos)

        for i, ftempl in enumerate(mode_ftempl_strs):
            if ftempl not in df_out.keys():
                continue
            redshiftTbl = df_out[ftempl]
            redChar = {
                'z_spec': np.array(redshiftTbl['z_spec']),
                'z_phot': np.array(redshiftTbl['z_phot']),
            }
            mask = (redChar['z_spec'] > zCharacteristic) & (~np.isnan(redChar['z_spec'])) & (~np.isnan(redChar['z_phot']) & (redChar['z_phot'] > 0))
            redChar['z_spec'] = redChar['z_spec'][mask]
            redChar['z_phot'] = redChar['z_phot'][mask]
            print(f"Have {len(redChar['z_spec'])} objects in {ftempl}")
            """deltaZ = (redChar['z_phot'] - redChar['z_spec'])/(1 + redChar['z_spec'])
            deltaZ_all = (redshiftTbl['z_phot'] - redshiftTbl['z_spec'])/(1 + redshiftTbl['z_spec'])
            outliers = (deltaZ > 0.15) | (deltaZ < -0.15)
            outliers_all = (deltaZ_all > 0.15) | (deltaZ_all < -0.15)
            bias = np.mean(deltaZ)
            scatter = np.std(deltaZ)"""
            deltaZ = deltaZ_s[mode][i]
            deltaZ_all = deltaZ_all_s[mode][i]
            outliers = outliers_s[mode][i]
            outliers_all = outliers_all_s[mode][i]
            bias = bias_s[mode][i]
            scatter = scatter_s[mode][i]
            eta = etas_s[mode][i]
            


            #outpath = outpaths.format(ftempl=ftempl, runTime=runTime)

            """#find light intensities
            ids = df_props['id'].values
            filtTab = Table.read(inpath, hdu=4)
            #get values ending with "CIRC1"
            pointIDs = filtTab['ID']
            filtTab = filtTab[[f for f in filtTab.colnames if f.endswith('CIRC1')]]
            #get dictionary version of table
            filtTab = filtTab.to_pandas()
            keys = filtTab.keys()
            fluxes = np.array([np.array(filtTab[key]) for key in keys]).T
            pointIntensities = np.sum(fluxes, axis=1)"""

            """#clear up bad ids and make simular sort
            for j,id in list(enumerate(ids))[::-1]:
                if id not in df_props['id'].values:
                    pointIntensities = np.delete(pointIntensities, j)
                    ids = np.delete(ids, j)
            #sort pointIntensities and ids to match df_props
            pointIntensities = pointIntensities[np.argsort(ids)]
            ids = np.sort(ids)
            antiSort = np.argsort(df_props['id'].values)
            sort = np.argsort(antiSort)
            pointIntensities = pointIntensities[sort]
            ids = ids[sort]"""
            #Remove sub zero and nan
            mask_cur = (redshiftTbl['z_spec'] > 0) & (~np.isnan(redshiftTbl['z_spec'])) & (~np.isnan(redshiftTbl['z_phot'])) & (redshiftTbl['z_phot'] > 0)
            x = np.array(redshiftTbl['z_spec'])
            y = np.array(redshiftTbl['z_phot'])
            #above = y > x + (1 + x) * 0.15 # outliers
            #below = y < x - (1 + x) * 0.15
            #outlier = above | below
            #mask_in = mask_cur & (~outlier)
            x = x[mask_cur]
            y = y[mask_cur]
            #mask_out = mask_cur & outlier
            #chi2_fit = redshiftTbl['z_phot_chi2'][mask_cur]/len(redshiftTbl.keys())#!new
            #avgchi2_fit = np.mean(chi2_fit)
            #medchi2_fit = np.median(chi2_fit)
            
            #avgchi2_fit_in = redshiftTbl['z_phot_chi2'][mask_in]/len(x)#!new
            #c = np.log10(pointIntensities[mask_cur])*((x/(1+x))**2)#!quilitative reshift scaling
            #c = np.ones(len(x))

            #calculate chi2 with linear regression of x and y
            chi2 = np.sum((y - x)**2/1)/len(x)#/len(x)#!dunno what to put for sigma

            xmin, xmax = np.c_[[x, y]].min(), np.c_[[x, y]].max()

            #axsMos[i].scatter(x, y, s=2.0, color='k')
            #hexbin with colorbar
            #axsMos[i].hexbin(x, y, gridsize=30, mincnt=1, edgecolors='none', bins='log')
            underCharc = (x < zCharacteristic)
            x_under, y_under = x[underCharc], y[underCharc]
            axsMos[i].hexbin(x_under, y_under, gridsize=30, mincnt=1, edgecolors='none', bins='log', cmap='viridis', alpha=0.2, extent=[xmin, xmax, xmin, xmax],zorder=200)
            axsMos[i].hexbin(redChar['z_spec'], redChar['z_phot'], gridsize=30, mincnt=1, edgecolors='none', bins='log', cmap='viridis', extent=[xmin, xmax, xmin, xmax])

            #plt.colorbar(axsMos[i].collections[0], ax=axsMos[i], label='log$_{10}$(N)')
            axsMos[i].plot([0, 100], [0, 100], c='r', ls='--', lw=0.5)
            axsMos[i].plot([zCharacteristic, zCharacteristic], [0, 14.5], c='r', ls='--', lw=0.5)
            
            #dict_stat = hmod.phot_spec_zs_stats(y, x)

            #annotate in top left
            annotAnchor = (0.05, 0.95)
            #axsMos[i].set_title(f'{ftempl}', fontsize=10)
            ax, ay = 2, 95        
            axsMos[i].annotate(f'{mode_ftempl_lbls[i]}', xy=annotAnchor, xycoords='axes fraction', fontsize=10, ha='left', va='top', xytext=(ax, ay), textcoords='axes points')

            #stats annotate
            sx, sy = 44, -55
            axsMos[i].annotate(f'$\eta$={round(eta,3):.3f}', xy=annotAnchor, xycoords='axes fraction', fontsize=8, ha='left', va='top', xytext=(ax+sx, ay+sy), textcoords='axes points')
            #etas = np.append(etas, eta**(-1))
            if mode_ftempl_lbls[i] == 'F45k' or mode_ftempl_lbls[i] == 'F60k': asterix = '*'
            else: asterix = ''
            #axsMos[i].annotate(f'$\eta_{{0+}}$={round(outliers_all.sum()/len(outliers_all),2):.2f}{asterix}', xy=annotAnchor, xycoords='axes fraction', fontsize=8, ha='left', va='top', xytext=(ax+sx, ay+sy-10), textcoords='axes points')
            axsMos[i].annotate(f'$\sigma_{{nmad}}$={round(np.median(np.abs(deltaZ[~outliers] - np.median(deltaZ[~outliers])))*1.4826,3):.3f}', xy=annotAnchor, xycoords='axes fraction', fontsize=8, ha='left', va='top', xytext=(ax+sx, ay+sy-10), textcoords='axes points')
            axsMos[i].annotate(f'$\Delta_{{bias}}$={round(bias,3):.3f}', xy=annotAnchor, xycoords='axes fraction', fontsize=8, ha='left', va='top', xytext=(ax+sx, ay+sy-20), textcoords='axes points')


            axsMos[i].set_xlim(zmin,zMax)
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