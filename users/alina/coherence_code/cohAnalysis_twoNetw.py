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
    fileName='CGplacebo_fix_nii_corrVals.pck'
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
    #titleName=condition+" coherence "
    #get3NetworkAvg(coherAvg_t, titleName, roiNames, numRuns)
    #titleName=condition+" correlation "
    #get3NetworkAvg(corrAvg_t, titleName, roiNames, numRuns)

    #Plot the data for 4 groups
    #Define the streams
    ventral=[1, 3, 6, 7, 8, 9, 10, 23, 25, 28, 29,  30, 31]
    dorsal=[14, 15, 16, 17, 18, 19, 20, 21, 35, 36, 37, 38, 39, 40, 41, 42]

    print 'Ventral rois: '+ str(roiNames[ventral])
    print 'Dorsal rois: ' + str(roiNames[dorsal])

    # Do Network analysis
    networkAvg= corrAll_t # Correlation (corrAll_t) or coherence (coherAll_t)
    analType='Correlation'

    # Get network averages
    ventCoher=getNetworkWithin(networkAvg, ventral)
    dorsCoher=getNetworkWithin(networkAvg, dorsal)

    #Average over last two dimensions....
    ventCoher_mean=stats.nanmean(ventCoher.reshape([numRuns,len(ventral)*len(ventral)]), axis=1)
    dorsCoher_mean=stats.nanmean(dorsCoher.reshape([numRuns,len(dorsal)*len(dorsal)]), axis=1)

    #Correlation means and STDs across all RUNs correlation/coherence values.
    allMeansWithin= (stats.nanmean(ventCoher_mean), stats.nanmean(dorsCoher_mean))
    allSTDWithin=(stats.nanstd(ventCoher_mean), stats.nanstd(dorsCoher_mean))

######
    # Get network btw
    #Early Visual
    EVbtwAllavg, EVbtwAllstd=getNetworkMeansBtw(networkAvg, ventral, [dorsal], numRuns)

    # Early Dorsal
    EDbtwAllavg, EDbtwAllstd=getNetworkMeansBtw(networkAvg, dorsal, [ventral], numRuns)
    allMeans=(np.insert(EVbtwAllavg, 0, allMeansWithin[0]), np.insert(EDbtwAllavg, 1, allMeansWithin[1]))
    allSTD=(np.insert(EVbtwAllstd, 0, allSTDWithin[0]), np.insert(EDbtwAllstd, 1, allSTDWithin[1]))

    # Make bar graph
    titleName=condition+ analType
    title= titleName+ ' by Network for ' +sub+ ' for '+ str(numRuns)+' runs'; labels=( 'Ventral', 'Dorsal')
    N=len(allMeansWithin)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.1       # the width of the bars
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, allMeans[0], width, color='r', yerr=allSTD[0])
    rects2 = ax.bar(ind+width*1, allMeans[1], width, color='b', yerr=allSTD[1])

    # add some labels
    ax.set_ylabel('Means')
    ax.set_title(title)
    ax.set_xticks(ind+width*2)
    ax.set_xticklabels( labels )
    ax.set_ylim( 0, 2.0 )
    ax.legend((rects1[0], rects2[0]),
        ('W/ Ventral', 'W/ Dorsal'))

    # Show final figure
    fig.show()

# Make a connection graph

#fig04 = drawgraph_channels(cohAll[sub], roiNames) #color_anchor=1
