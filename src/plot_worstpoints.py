import numpy as np
import matplotlib.pyplot as plt

def plot(df_out, ftempl_strs, runTime):
    
    redshiftPoints = np.sum([df_out[ftempl_strs[i]]['z_phot'] for i in range(len(ftempl_strs))],axis=0) / len(ftempl_strs)
    redshiftPointsSpec = np.sum([df_out[ftempl_strs[i]]['z_spec'] for i in range(len(ftempl_strs))],axis=0) / len(ftempl_strs)
    redshiftPointsErr = np.std([df_out[ftempl_strs[i]]['z_phot'] for i in range(len(ftempl_strs))],axis=0)
    redshiftPointsDelta = np.abs(redshiftPoints - redshiftPointsSpec)
    redshiftPointsIDs = df_out[ftempl_strs[0]]['ID']

    #sort by delta
    sortIdx = np.argsort(redshiftPointsDelta/redshiftPointsErr)[::-1]
    redshiftPoints = redshiftPoints[sortIdx]
    redshiftPointsSpec = redshiftPointsSpec[sortIdx]
    redshiftPointsErr = redshiftPointsErr[sortIdx]
    redshiftPointsDelta = redshiftPointsDelta[sortIdx]
    redshiftPointsIDs = redshiftPointsIDs[sortIdx]

    #print first 10
    idsOfInterest = []
    print("Worst 10 points:")
    fig, ax = plt.subplots(figsize=(10,2))
    c=0
    maxZSpec = 0
    maxZPhot = 0
    yTicks = []
    for i in range(len(redshiftPoints)):
        delta = redshiftPointsDelta[i]
        sigma_phot = redshiftPointsErr[i]
        z_spec = redshiftPointsSpec[i]
        z_phot = redshiftPoints[i]
        id = redshiftPointsIDs[i]
        if delta-(sigma_phot*1.3) < 0.0: continue
        if z_spec < 0.0: continue
        if z_phot < 0.0: continue
        if delta < 1: continue
        #if z_spec < 4.0: continue
        print(f"z_spec: {z_spec:.3f}, z_phot: {z_phot:.3f}, delta: {delta:.3f}, id: {id}, sigma_phot: {sigma_phot:.3f}")
        
        label = None
        if c == 0: label = "z_phot - all templates"
        plt.errorbar(z_phot, -c, xerr=sigma_phot, fmt='x', color='k', capsize=2, label=label)
        if z_phot > maxZPhot: maxZPhot = z_phot
        label = None
        if c == 0: label = "z_spec"
        plt.plot([z_spec,z_spec],[-c-0.5,-c+0.5],color='r',linestyle='--',label=label,linewidth=2)
        if z_spec > maxZSpec: maxZSpec = z_spec
        plt.plot([0,100],[-c,-c],color='k',linestyle='-',linewidth=1,alpha=0.1)

        yTicks.append(id)

        if sigma_phot < 0.5: 
            idsOfInterest.append(id)
            print("Added to list of interest")
        c+=1

    #set y ticks to c locations
    ax.set_yticks(-np.arange(len(yTicks)))
    #remove numbers on yticks
    ax.set_yticklabels([])
    #annotate with IDs on yticks
    for i,id in enumerate(yTicks):
        ax.annotate(str(id) + "  ", xy=(0, -i+0.7), xycoords='data', fontsize=6, ha='right', va='top',rotation=20)

    plt.xlim(0,max(maxZSpec+0.5,maxZPhot+1))
    plt.xlabel("Redshift (z)")
    plt.ylabel("Object")
    ax.yaxis.set_label_coords(-0.05,0.5)

    #set upper margin to 0.01
    #ax.set_ylim(-c-0.5,2)
    plt.legend(loc='lower right')

    #savefig
    plt.savefig(f"figures/forpaper/consistent_outliers_{runTime}.png",dpi=300,bbox_inches='tight')

    plt.show()