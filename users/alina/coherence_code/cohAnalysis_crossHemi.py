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

def getNetworkMeansBtw(networkAvg, netw1, other_netw, numRuns):
    allMeans=np.zeros(len(other_netw))
    allSTDs=np.zeros(len(other_netw))
    i=0
    for net in other_netw:
            netw1btwnetw_mean=getNetworkBtw(networkAvg, netw1, net, numRuns)
            allMeans[i]=np.mean(netw1btwnetw_mean)
            allSTDs[i]=np.std(netw1btwnetw_mean)
            i=i+1

    return allMeans, allSTDs

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

    networkAvg= corrAll_t# Correlation (corrAll_t) or coherence (coherAll_t)
    analType='Correlation'

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
    rEVbtwAllavg, rEVbtwAllstd=getNetworkMeansBtw(networkAvg, r_earlyVent, [r_earlyDors, r_parietal, r_objSel, l_earlyVent, l_earlyDors, l_parietal, l_objSel], numRuns)
    lEVbtwAllavg, lEVbtwAllstd=getNetworkMeansBtw(networkAvg, l_earlyVent, [r_earlyVent, r_earlyDors,  r_parietal, r_objSel, l_earlyDors, l_parietal, l_objSel], numRuns)

    # Early Dorsal
    rEDbtwAllavg, rEDbtwAllstd=getNetworkMeansBtw(networkAvg, r_earlyDors, [r_earlyVent, r_parietal, r_objSel, l_earlyVent, l_earlyDors, l_parietal, l_objSel], numRuns)
    lEDbtwAllavg, lEDbtwAllstd=getNetworkMeansBtw(networkAvg, l_earlyDors, [r_earlyVent, r_earlyDors, r_parietal, r_objSel, l_earlyVent, l_parietal, l_objSel], numRuns)

    # Parietal
    rParbtwAllavg, rParbtwAllstd=getNetworkMeansBtw(networkAvg, r_parietal, [r_earlyVent, r_earlyDors, r_objSel, l_earlyVent, l_earlyDors, l_objSel, l_parietal], numRuns)
    lParbtwAllavg, lParbtwAllstd=getNetworkMeansBtw(networkAvg, l_parietal, [r_earlyVent, r_earlyDors, r_parietal, r_objSel, l_earlyVent, l_earlyDors, l_objSel ], numRuns)

    # Object Selective
    rObjSelbtwAllavg, rObjSelbtwAllstd=getNetworkMeansBtw(networkAvg, r_objSel, [r_earlyVent, r_earlyDors, r_parietal, l_earlyVent, l_earlyDors, l_parietal, l_objSel], numRuns)
    lObjSelbtwAllavg, lObjSelbtwAllstd=getNetworkMeansBtw(networkAvg, l_objSel, [r_earlyVent, r_earlyDors, r_parietal, r_objSel, l_earlyVent, l_earlyDors, l_parietal], numRuns)


    allMeans=(np.insert(rEVbtwAllavg, 0, allMeansWithin[0]), np.insert(rEDbtwAllavg, 1, allMeansWithin[1]), np.insert(rParbtwAllavg, 2, allMeansWithin[2]),
        np.insert(rObjSelbtwAllavg, 3, allMeansWithin[3]), np.insert(lEVbtwAllavg, 4, allMeansWithin[4]), np.insert(lEDbtwAllavg, 5, allMeansWithin[5]),
        np.insert(lParbtwAllavg, 6, allMeansWithin[6]), np.insert(lObjSelbtwAllavg, 7, allMeansWithin[7]))
    allSTD=(np.insert(rEVbtwAllstd, 0, allSTDWithin[0]), np.insert(rEDbtwAllstd, 1, allSTDWithin[1]), np.insert(rParbtwAllstd, 2, allSTDWithin[2]),
        np.insert(rObjSelbtwAllstd, 3, allSTDWithin[3]), np.insert(lEVbtwAllstd, 4, allSTDWithin[4]), np.insert(lEDbtwAllstd, 5, allSTDWithin[5]),
        np.insert(lParbtwAllstd, 6, allSTDWithin[6]), np.insert(lObjSelbtwAllstd, 7, allSTDWithin[7]))


    titleName=condition+ analType
    # Make bar graph
    title= titleName+ ' by Network for ' +sub+ ' for '+ str(numRuns)+' runs'; labels=( 'rEarly Ventral', 'rEarly Dorsal', 'rParietal', 'rObject Selective', 'lEarly Ventral', 'lEarly Dorsal', 'lParietal', 'lObject Selective')

    N=len(allMeansWithin)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.1       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, allMeans[0], width, color='r', yerr=allSTD[0])
    rects2 = ax.bar(ind+width*1, allMeans[1], width, color='y', yerr=allSTD[1])
    rects3 = ax.bar(ind+width*2, allMeans[2], width, color='g', yerr=allSTD[2])
    rects4 = ax.bar(ind+width*3, allMeans[3], width, color='b', yerr=allSTD[3])
    rects5 = ax.bar(ind+width*4, allMeans[4], width, color='r', yerr=allSTD[4])
    rects6 = ax.bar(ind+width*5, allMeans[5], width, color='y', yerr=allSTD[5])
    rects7 = ax.bar(ind+width*6, allMeans[6], width, color='g', yerr=allSTD[6])
    rects8 = ax.bar(ind+width*7, allMeans[7], width, color='b', yerr=allSTD[7])



    # add some labels
    ax.set_ylabel('Means')
    ax.set_title(title)
    ax.set_xticks(ind+width*2)
    ax.set_xticklabels( labels )
    ax.set_ylim( 0, 2.0 )
    ax.legend((rects1[0], rects2[0], rects3[0], rects4[0], rects5[0], rects6[0], rects7[0], rects8[0]),
        ('W/ rEarly Ventral', 'W/ rEarly Dorsal', 'W/ rParietal', 'W/ rObject Sel.', 'W/ lEarly Ventral', 'W/ lEarly Dorsal', 'W/ lParietal', 'W/ lObject Sel.'))

    # Show final figure
    fig.show()

# Make a connection graph

#fig04 = drawgraph_channels(cohAll[sub], roiNames) #color_anchor=1
