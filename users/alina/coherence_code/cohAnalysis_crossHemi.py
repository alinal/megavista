#Changed on 1/30
import numpy as np
import scipy.stats as stats
import pickle
import datetime

# Import network analysis functions
import networkAnal
reload(networkAnal)
from networkAnal import get3NetworkAvg, getNetworkWithin, getNetworkBtw


from matplotlib import mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from mpl_toolkits.axes_grid import make_axes_locatable

import vista_utils as tsv # Get it at: https://github.com/arokem/vista_utils
from nitime.fmri.io import time_series_from_file as load_nii
import nitime.timeseries as ts
import nitime.viz as viz

from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr

#Some visualization functions require networkx. Import that if possible:
try:
    import networkx as nx
#If not, throw an error and get on with business:
except ImportError:
    e_s = "Networkx is not available. Some visualization tools might not work"
    e_s += "\n To download networkx: http://networkx.lanl.gov/"
    print e_s
    class NetworkxNotInstalled(object):
        def __getattribute__(self,x):
            raise ImportError(e_s)
            nx = NetworkxNotInstalled()

def makeBarPlots(allMeansWithin, allSTDWithin, allMeansBtw, allSTDBtw, title, labels):
    N=len(allMeansWithin)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, allMeansWithin, width, color='r', yerr=allSTDWithin)

    rects2 = ax.bar(ind+width, allMeansBtw, width, color='y', yerr=allSTDBtw)

    # add some
    ax.set_ylabel('Means')
    ax.set_title(title)
    ax.set_xticks(ind+width)
    ax.set_xticklabels( labels )
    ax.set_ylim( 0, 1.5 )

    ax.legend( (rects1[0], rects2[0]), ('Within', 'Between') )

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height*100),
            ha='center', va='bottom')

        autolabel(rects1)
        autolabel(rects2)

        plt.show()


if __name__ == "__main__":
    # Close any opened plots
    plt.close('all')

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'
    fileName='CG&CHT&DCA&SSvizROIsOrderfix_normalize_bothplacebo1runs_2013-08-02.pck'
    #fileName='CG&CHT&DCAallROIsOrderFix_normalizeplacebo1runs_2012-02-08.pck'
    #fileName='CG&CHT&DCAallROIsOrderFix_normalizedonepazil1runs_2013-05-29.pck'
    #fileName='CG&CHT&DCAallROIsOrderLeft_normalizeplacebo1runs_2013-05-29.pck'
    condition='NormalizedFixationRun'
    loadFile=fmri_path+'Results/' +fileName

    figSize=[10., 10.]
    # Load the data
    print 'Loading subject coherence and correlation dictionaries.'
    file=open(loadFile, 'r')
    # First file loaded is coherence
    cohAll=pickle.load(file)
    # Second file loaded is correlation
    corrAll=pickle.load(file)
    roiNames=pickle.load(file)
    file.close()

for sub in cohAll:
    print sub
    numRuns=cohAll[sub].shape[0]

    #Fisher transform the data (maybe fisher transform earlier)
    coherAll_t = np.arctanh(cohAll[sub][:])
    corrAll_t=np.arctanh(corrAll[sub][:])

    #replace all inf in fisher transform with nan
    coherAll_flat=coherAll_t.flatten()
    corrAll_flat=corrAll_t.flatten()
    ind = np.where(coherAll_flat == np.Infinity)
    coherAll_flat[ind] = np.nan
    ind = np.where( corrAll_flat == np.Infinity)
    corrAll_flat[ind] = np.nan

    # Reshape back into matrix of runs x rois x rois
    coherAll_t=coherAll_flat.reshape(numRuns, roiNames.size, roiNames.size)
    corrAll_t=corrAll_flat.reshape(numRuns, roiNames.size, roiNames.size)

    # Average over runs (the first dimension)
    coherAvg_t=np.mean(coherAll_t, 0)
    coherSTD=np.std(coherAll_t, 0)
    corrAvg_t=np.mean(corrAll_t,0)
    corrSTD=np.std(corrAll_t, 0)

    # Plot graph of coherence and correlation values
    fig1 = drawmatrix_channels(coherAvg_t, roiNames, size=coherAvg_t.shape, color_anchor=0, title='Average ' +condition+  ' Coherence Results over ' +str(numRuns) + ' runs for ' + sub)
    fig2=drawmatrix_channels(coherSTD, roiNames, size=coherSTD.shape, color_anchor=0, title='Average ' +condition+ ' Coherence STD over ' +str(numRuns) + ' runs for ' + sub)
    fig3=drawmatrix_channels(corrAvg_t, roiNames, size=corrAvg_t.shape, color_anchor=0, title='Average ' +condition+ ' Correlation Results over ' +str(numRuns) + ' runs for ' + sub)
    fig4=drawmatrix_channels(corrSTD, roiNames, size=corrSTD.shape, color_anchor=0, title='Average ' +condition+ ' Correlation STD over ' +str(numRuns) + ' runs for ' + sub)
    plt.show()

    #Plot data for 3 streams (btw for all)
    titleName=condition+" coherence "
    #get3NetworkAvg(coherAvg_t, titleName, roiNames, numRuns)
    titleName=condition+" correlation "
    #get3NetworkAvg(corrAvg_t, titleName, roiNames, numRuns)

    #Plot the data for 4 groups
    #Define the streams
    r_earlyVent=[1, 2, 3]
    r_earlyDors=[7, 8, 9]
    r_parietal=[10, 11, 12]
    r_objSel=[4, 5, 6]
    l_earlyVent=[14, 15, 16]
    l_earlyDors=[19, 20, 21]
    l_parietal=[ 22, 23, 24]
    l_objSel=[ 17, 18]

    print 'Early Ventral rois: '+ str(roiNames[r_earlyVent]) + str(roiNames[l_earlyVent])
    print 'Early Dorsal rois: ' + str(roiNames[r_earlyDors])+ str(roiNames[l_earlyDors])
    print 'Parietal rois: '+ str(roiNames[r_parietal]) +str(roiNames[l_parietal])
    print 'Object rois: '+ str(roiNames[r_objSel])+ str(roiNames[l_objSel])

    networkAvg= coherAll_t# Correlation (corrAll_t) or coherence (coherAll_t)
    analType='Coherence'

    # Get network averages
    r_earlyVentCoher=getNetworkWithin(networkAvg, r_earlyVent)
    r_earlyDorsCoher=getNetworkWithin(networkAvg, r_earlyDors)
    r_parietalCoher=getNetworkWithin(networkAvg, r_parietal)
    r_objSelCoher=getNetworkWithin(networkAvg, r_objSel)
    l_earlyVentCoher=getNetworkWithin(networkAvg, l_earlyVent)
    l_earlyDorsCoher=getNetworkWithin(networkAvg, l_earlyDors)
    l_parietalCoher=getNetworkWithin(networkAvg, l_parietal)
    l_objSelCoher=getNetworkWithin(networkAvg, l_objSel)

    #Average over last two dimensions....
    r_earlyVentCoher_mean=stats.nanmean(r_earlyVentCoher.reshape([numRuns,len(r_earlyVent)*len(r_earlyVent)]), axis=1)
    r_earlyDorsCoher_mean=stats.nanmean(r_earlyDorsCoher.reshape([numRuns,len(r_earlyDors)*len(r_earlyDors)]), axis=1)
    r_parietalCoher_mean=stats.nanmean(r_parietalCoher.reshape([numRuns,len(r_parietal)*len(r_parietal)]), axis=1)
    r_objSelCoher_mean=stats.nanmean(r_objSelCoher.reshape([numRuns,len(r_objSel)*len(r_objSel)]), axis=1)
    l_earlyVentCoher_mean=stats.nanmean(l_earlyVentCoher.reshape([numRuns,len(l_earlyVent)*len(l_earlyVent)]), axis=1)
    l_earlyDorsCoher_mean=stats.nanmean(l_earlyDorsCoher.reshape([numRuns,len(l_earlyDors)*len(l_earlyDors)]), axis=1)
    l_parietalCoher_mean=stats.nanmean(l_parietalCoher.reshape([numRuns,len(l_parietal)*len(l_parietal)]), axis=1)
    l_objSelCoher_mean=stats.nanmean(l_objSelCoher.reshape([numRuns,len(l_objSel)*len(l_objSel)]), axis=1)

    #Correlation means and STDs across all RUNs correlation/coherence values.
    allMeansWithin= (stats.nanmean(r_earlyVentCoher_mean), stats.nanmean(r_earlyDorsCoher_mean), stats.nanmean(r_parietalCoher_mean), stats.nanmean(r_objSelCoher_mean),
        stats.nanmean(l_earlyVentCoher_mean), stats.nanmean(l_earlyDorsCoher_mean), stats.nanmean(l_parietalCoher_mean), stats.nanmean(l_objSelCoher_mean) )
    allSTDWithin=(stats.nanstd(r_earlyVentCoher_mean), stats.nanstd(r_earlyDorsCoher_mean), stats.nanstd(r_parietalCoher_mean), stats.nanstd(r_objSelCoher_mean),
        stats.nanstd(l_earlyVentCoher_mean), stats.nanstd(l_earlyDorsCoher_mean), stats.nanstd(l_parietalCoher_mean), stats.nanstd(l_objSelCoher_mean))

######
    # Get network btw
    #Early Visual
    EVbtwED_mean=getNetworkBtw(networkAvg, earlyVent, earlyDors, numRuns); EVbtwEDavg=np.mean(EVbtwED_mean); EVbtwEDstd=np.std(EVbtwED_mean)
    EVbtwPar_mean=getNetworkBtw(networkAvg, earlyVent, parietal, numRuns); EVbtwParavg=np.mean(EVbtwPar_mean); EVbtwParstd=np.std(EVbtwPar_mean)
    EVbtwObjSel_mean=getNetworkBtw(networkAvg, earlyVent, objSel, numRuns); EVbtwObjSelavg=np.mean(EVbtwObjSel_mean); EVbtwObjSelstd=np.std(EVbtwObjSel_mean)

    # Early Dorsal
    EDbtwEV_mean=getNetworkBtw(networkAvg, earlyDors, earlyVent, numRuns); EDbtwEVavg=np.mean(EDbtwEV_mean); EDbtwEVstd=np.std(EDbtwEV_mean)
    EDbtwPar_mean=getNetworkBtw(networkAvg, earlyDors, parietal, numRuns); EDbtwParavg=np.mean(EDbtwPar_mean); EDbtwParstd=np.std(EDbtwPar_mean)
    EDbtwObjSel_mean=getNetworkBtw(networkAvg, earlyDors, objSel, numRuns); EDbtwObjSelavg=np.mean(EDbtwObjSel_mean); EDbtwObjSelstd=np.std(EDbtwObjSel_mean)

    # Parietal
    ParbtwEV_mean=getNetworkBtw(networkAvg, parietal, earlyVent, numRuns); ParbtwEVavg=np.mean(ParbtwEV_mean); ParbtwEVstd=np.std(ParbtwEV_mean)
    ParbtwED_mean=getNetworkBtw(networkAvg, parietal, earlyDors, numRuns); ParbtwEDavg=np.mean(ParbtwED_mean); ParbtwEDstd=np.std(ParbtwED_mean)
    ParbtwObjSel_mean=getNetworkBtw(networkAvg, parietal, objSel, numRuns); ParbtwObjSelavg=np.mean(ParbtwObjSel_mean); ParbtwObjSelstd=np.std(ParbtwObjSel_mean)

    # Object Selective
    ObjSelbtwEV_mean=getNetworkBtw(networkAvg, objSel, earlyVent, numRuns); ObjSelbtwEVavg=np.mean(ObjSelbtwEV_mean); ObjSelbtwEVstd=np.std(ObjSelbtwEV_mean)
    ObjSelbtwED_mean=getNetworkBtw(networkAvg, objSel, earlyDors, numRuns); ObjSelbtwEDavg=np.mean(ObjSelbtwED_mean); ObjSelbtwEDstd=np.std(ObjSelbtwED_mean)
    ObjSelbtwPar_mean=getNetworkBtw(networkAvg, objSel, parietal, numRuns); ObjSelbtwParavg=np.mean(ObjSelbtwPar_mean); ObjSelbtwParstd=np.std(ObjSelbtwPar_mean)


    allMeans=([allMeansWithin[0], EVbtwEDavg, EVbtwParavg, EVbtwObjSelavg], [EDbtwEVavg, allMeansWithin[1], EDbtwParavg, EDbtwObjSelavg],
        [ParbtwEVavg, ParbtwEDavg, allMeansWithin[2], ParbtwObjSelavg], [ObjSelbtwEVavg, ObjSelbtwEDavg, ObjSelbtwParavg, allMeansWithin[3]])
    allSTD=([allSTDWithin[0], EVbtwEDstd, EVbtwParstd, EVbtwObjSelstd], [EDbtwEVstd, allSTDWithin[1], EDbtwParstd, EDbtwObjSelstd], [ParbtwEVstd, ParbtwEDstd, allSTDWithin[2], ParbtwObjSelstd],
        [ObjSelbtwEVstd, ObjSelbtwEDstd, ObjSelbtwParstd, allSTDWithin[3]])

    titleName=condition+ analType
    # Make bar graph
    title= titleName+ ' by Network for ' +sub+ ' for '+ str(numRuns)+' runs'; labels=( 'Early Ventral', 'Early Dorsal', 'Parietal', 'Object Selective')

    N=len(allMeansWithin)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, allMeans[0], width, color='r', yerr=allSTD[0])
    rects2 = ax.bar(ind+width*1, allMeans[1], width, color='y', yerr=allSTD[1])
    rects3 = ax.bar(ind+width*2, allMeans[2], width, color='g', yerr=allSTD[2])
    rects4 = ax.bar(ind+width*3, allMeans[3], width, color='b', yerr=allSTD[3])

    # add some labels
    ax.set_ylabel('Means')
    ax.set_title(title)
    ax.set_xticks(ind+width*2)
    ax.set_xticklabels( labels )
    ax.set_ylim( 0, 2.0 )
    ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('W/ Early Ventral', 'W/ Early Dorsal', 'W/ Parietal', 'W/ Object Sel.'))

    # Show final figure
    fig.show()

# Make a connection graph

#fig04 = drawgraph_channels(cohAll[sub], roiNames) #color_anchor=1
