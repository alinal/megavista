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
reload(viz)

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

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
    #plt.close('all')
    # Plot Matrix
    plotMat=0

    base_path = '/Users/Alina/Desktop/' # Change this to your path
    fmri_path = base_path + 'fmri/'
    condition='Fixation'
    fileName='CGplacebo_fix_nii_43ROIts_corrVals_wGM_hierarch_22reg.pck'
    #fileName='CG&CHT&DCAallROIsOrderFix_normalizeplacebo1runs_2012-02-08.pck'
    #fileName='CG&CHT&DCAallROIsOrderFix_normalizedonepazil1runs_2013-05-29.pck'
    #fileName='CG&CHT&DCAallROIsOrderLeft_normalizeplacebo1runs_2013-05-29.pck'

    #loadFile=fmri_path+'Results/correlation/' +fileName
    loadFile=base_path+fileName

    figSize=[10., 10.]
    # Load the data
    print 'Loading subject coherence and correlation dictionaries from %s.' % fileName
    file=open(loadFile, 'r')
    # First file loaded is coherence
    cohAll=pickle.load(file)
    # Second file loaded is correlation
    corrAll=pickle.load(file)
    roiNames=pickle.load(file)
    subs=pickle.load(file)
    hierarch=pickle.load(file)
    file.close()

for sub in cohAll:
    print sub
    numRuns=cohAll[sub].shape[0]

    #  Do hierarchical clustering of all ROIs
    hierarch_t=np.mean(hierarch[sub],0)
    Z=linkage(hierarch_t, 'single')
    plt.figure()
    dendrogram(Z, color_threshold=0, labels=roiNames)

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

    if plotMat:
        # Plot graph of coherence and correlation values
        fig1 = drawmatrix_channels(coherAvg_t, roiNames, size=coherAvg_t.shape, color_anchor=0) # title='Average ' +condition+  ' Coherence Results over ' +str(numRuns) + ' runs for ' + sub)
        fig2=drawmatrix_channels(coherSTD, roiNames, size=coherSTD.shape, color_anchor=0, title='Average ' +condition+ ' Coherence STD over ' +str(numRuns) + ' runs for ' + sub)
        fig3=drawmatrix_channels(corrAvg_t, roiNames, size=corrAvg_t.shape, color_anchor=0) # title='Average ' +condition+ ' Correlation Results over ' +str(numRuns) + ' runs for ' + sub)
        fig4=drawmatrix_channels(corrSTD, roiNames, size=corrSTD.shape, color_anchor=0, title='Average ' +condition+ ' Correlation STD over ' +str(numRuns) + ' runs for ' + sub)
        plt.show()

    #Define the streams
    ventralRH=['R_V2V_0.25', 'R_V2D_0.25', 'R_V3V_0.25', 'R_V3D_0.25', 'R_V4_0.25', 'r_IOG_p3_0.25', 'r_LOf_p3_0.25', 'r_pFus_p3_0.25', 'r_mFus_p3_0.25', 'r_PPA_p4_0.25']
    ventralLH=['L_V2V_0.25', 'L_V2D_0.25', 'L_V3V_0.25', 'L_V3D_0.25', 'L_V4_0.25', 'l_IOG_p3_0.25', 'l_LOf_p3_0.25', 'l_pFus_p3_0.25', 'l_mFus_p3_0.25', 'l_PPA_p4_0.25']
    dorsalRH=['R_V3A_0.25', 'R_MT_al_.5_0.25', 'R_IPS0_0.25', 'R_IPS1_0.25', 'R_IPS2_0.25', 'R_IPS3_0.25', 'R_IPS4_0.25', 'R_IPS5_0.25']
    dorsalLH=[ 'L_V3A_0.25', 'L_MT_al_.5_0.25', 'L_IPS0_0.25', 'L_IPS1_0.25', 'L_IPS2_0.25', 'L_IPS3_0.25', 'L_IPS4_0.25', 'L_IPS5_0.25']

    ventralRHIndx=np.where(np.in1d(roiNames, ventralRH))[0]
    ventralLHIndx=np.where(np.in1d(roiNames, ventralLH))[0]
    dorsalRHIndx=np.where(np.in1d(roiNames, dorsalRH))[0]
    dorsalLHIndx=np.where(np.in1d(roiNames, dorsalLH))[0]

    print 'Ventral RH rois: '+ str(roiNames[ventralRHIndx])
    print 'Ventral LH rois: '+ str(roiNames[ventralLHIndx])
    print 'Dorsal RH rois: ' + str(roiNames[dorsalRHIndx])
    print 'Dorsal LH rois: ' + str(roiNames[dorsalLHIndx])

    # Do Network analysis
    networkAvg= corrAll_t # Correlation (corrAll_t) or coherence (coherAll_t)
    analType='Correlation'

    # Get network averages
    ventCoherRH=getNetworkWithin(networkAvg, ventralRHIndx)
    ventCoherLH=getNetworkWithin(networkAvg, ventralLHIndx)
    dorsCoherRH=getNetworkWithin(networkAvg, dorsalRHIndx)
    dorsCoherLH=getNetworkWithin(networkAvg, dorsalLHIndx)

    #Average over last two dimensions....
    ventCoherRH_mean=stats.nanmean(ventCoherRH.reshape([numRuns,len(ventralRHIndx)*len(ventralRHIndx)]), axis=1)
    dorsCoherRH_mean=stats.nanmean(dorsCoherRH.reshape([numRuns,len(dorsalRHIndx)*len(dorsalRHIndx)]), axis=1)
    ventCoherLH_mean=stats.nanmean(ventCoherLH.reshape([numRuns,len(ventralLHIndx)*len(ventralLHIndx)]), axis=1)
    dorsCoherLH_mean=stats.nanmean(dorsCoherLH.reshape([numRuns,len(dorsalLHIndx)*len(dorsalLHIndx)]), axis=1)

    #Correlation means and STDs across all RUNs correlation/coherence values.
    allMeansWithin= (stats.nanmean(ventCoherRH_mean), stats.nanmean(ventCoherLH_mean), stats.nanmean(dorsCoherRH_mean), stats.nanmean(dorsCoherLH_mean))
    allSTDWithin=(stats.nanstd(ventCoherRH_mean), stats.nanstd(ventCoherLH_mean), stats.nanstd(dorsCoherRH_mean), stats.nanstd(dorsCoherLH_mean))

######
    # Get network btw
    #Early Visual
    rhEVbtwAllavg, rhEVbtwAllstd=getNetworkMeansBtw(networkAvg, ventralRHIndx, [ventralLHIndx, dorsalRHIndx, dorsalLHIndx], numRuns)
    lhEVbtwAllavg, lhEVbtwAllstd=getNetworkMeansBtw(networkAvg, ventralLHIndx, [ventralRHIndx, dorsalRHIndx, dorsalLHIndx], numRuns)

    # Early Dorsal
    rhEDbtwAllavg, rhEDbtwAllstd=getNetworkMeansBtw(networkAvg, dorsalRHIndx, [ventralRHIndx, ventralLHIndx, dorsalLHIndx], numRuns)
    lhEDbtwAllavg, lhEDbtwAllstd=getNetworkMeansBtw(networkAvg, dorsalLHIndx, [ventralRHIndx, ventralLHIndx, dorsalRHIndx], numRuns)

    allMeans=(np.insert(rhEVbtwAllavg, 0, allMeansWithin[0]),np.insert(lhEVbtwAllavg, 1, allMeansWithin[1]), np.insert(rhEDbtwAllavg, 2, allMeansWithin[2]), np.insert(lhEDbtwAllavg, 3, allMeansWithin[3]))
    allSTD=(np.insert(rhEVbtwAllstd, 0, allSTDWithin[0]),np.insert(lhEVbtwAllstd, 1, allSTDWithin[1]), np.insert(rhEDbtwAllstd, 2, allSTDWithin[2]), np.insert(lhEDbtwAllstd, 3, allSTDWithin[3]))
    allSTD=allSTD/np.sqrt(numRuns)

    # Make bar graph
    titleName=condition+ ' ' + analType
    title= titleName+ ' by Network for ' +sub+ ' for '+ str(numRuns)+' runs'; labels=( 'rhVentral', 'lhVentral', 'rhDorsal', 'lhDorsal')
    N=len(allMeansWithin)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.1       # the width of the bars
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, allMeans[0], width, color='b', yerr=allSTD[0])
    rects2 = ax.bar(ind+width*1, allMeans[1], width, color='r', yerr=allSTD[1])
    rects3 = ax.bar(ind+width*2, allMeans[2], width, color='g', yerr=allSTD[2])
    rects4 = ax.bar(ind+width*3, allMeans[3], width, color='m', yerr=allSTD[3])

    # Print results
    print 'All means rh ventral: %s' % allMeans[0]
    print 'All means lh ventral: %s' % allMeans[1]
    print 'All means rh dorsal: %s' % allMeans[2]
    print 'All means lh dorsal: %s' % allMeans[3]


    # add some labels
    ax.set_ylabel('Means')
    ax.set_title(title)
    ax.set_xticks(ind+width*2)
    ax.set_xticklabels( labels )
    ax.set_ylim( 0, 1.0 )
    ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
        ('W/ rh Ventral',  'W/ lh Ventral', 'W/ rh Dorsal', 'W/ lh Dorsal'))

    # Show final figure
    fig.show()
    1/0

# Make a connection graph
allIndx=np.concatenate([ventralRHIndx, dorsalRHIndx])
rhROIs=coherAvg_t[allIndx,:][:,allIndx]
print roiNames[allIndx]

#how to do the graph:
do_tstat = False
thresh = .75 #abs thresh for graph (either t or fisher-r, depending on do_tstat), 0 = no thresh
    
#make into a t-stat map if called upon
if do_tstat:
    rhROIs = cmat[net]/cmat_se[net]
    nweight = 2
else:
    nweight = 25

plt.figure(figsize=(8,8))

# Zero out values below diagonal
rhROIs=np.triu(rhROIs, k=1)

# Make a graph object
G1 = nx.Graph(weighted = True)
G1 = nx.from_numpy_matrix(rhROIs,G1)

# Set up labels and partitions
nnod=np.shape(rhROIs)[1]
nod_labels=roiNames[allIndx]

pos=nx.circular_layout(G1)

# Add nodes to the plot
nx.draw_networkx_nodes(G1,pos,alpha=0.5,node_color='w')

#draw positive edges
evals_pos = np.array([d['weight'] for (u,v,d) in G1.edges(data=True) if d['weight']>thresh])
e_pos = [(u,v) for (u,v,d) in G1.edges(data=True) if d['weight']>thresh]

nx.draw_networkx_edges(G1,pos,edgelist=e_pos,width=evals_pos*nweight,alpha=1,edge_color='r')
                    
#draw negative edges
evals_neg = np.array([d['weight'] for (u,v,d) in G1.edges(data=True) if d['weight']<-1*thresh])
e_neg = [(u,v) for (u,v,d) in G1.edges(data=True) if d['weight']<-1*thresh]
        
nx.draw_networkx_edges(G1,pos,edgelist=e_neg,width=evals_neg*-1*nweight,alpha=1,edge_color='b')

#draw labels
nx.draw_networkx_labels(G1,pos,nod_labels,font_size=8,font_weight='bold')

#fig04 = drawgraph_channels(cohAll[sub], roiNames) #color_anchor=1
