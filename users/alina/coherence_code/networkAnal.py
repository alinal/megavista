from nitime.utils import triu_indices
import numpy as np



def get3NetworkAvg(data_t, titleName, roiNames, numRuns):
    #Define the streams
    #Ventral=[1, 3, 11, 12, 13, 14]
    #Dorsal=[2, 4, 5, 6, 7, 8, 9, 10]
    #Lateral=[0, 1, 2, 3, 4]
    Lateral=[0,1,2,8,9]
    Dorsal=[8,9,10, 11, 12, 13, 14, 15]
    Ventral=[1,2, 3, 4, 5, 6]

    print 'Ventral rois: '+ str(roiNames[Ventral])
    print 'Dorsal rois: ' + str(roiNames[Dorsal])
    print 'Early Visual rois: '+ str(roiNames[Lateral])

    # Get network averages
    lateralCoher=getNetworkWithin(data_t, Lateral)
    dorsalCoher=getNetworkWithin(data_t, Ventral)
    ventralCoher=getNetworkWithin(data_t, Dorsal)
    #allMeansWithin=(stats.nanmean(lateralCoher.flat), stats.nanmean(dorsalCoher.flat), stats.nanmean(ventralCoher.flat))
    #allSTDWithin=(stats.nanstd(lateralCoher.flat), stats.nanstd(dorsalCoher.flat), stats.nanstd(ventralCoher.flat))
    allMeansWithin= (stats.nanmean(dorsalCoher.flat), stats.nanmean(ventralCoher.flat))
    allSTDWithin=( stats.nanstd(dorsalCoher.flat), stats.nanstd(ventralCoher.flat))

    latBtwCoher=getNetworkBtw(data_t, Lateral, Ventral+Dorsal)
    dorsBtwCoher=getNetworkBtw(data_t, Dorsal, Ventral)
    ventBtwCoher=getNetworkBtw(data_t, Ventral, Dorsal)

    #allMeansBtw=(stats.nanmean(latBtwCoher), stats.nanmean(dorsBtwCoher), stats.nanmean(ventBtwCoher))
    #allSTDBtw=(stats.nanstd(latBtwCoher), stats.nanstd(dorsBtwCoher), stats.nanstd(ventBtwCoher))
    # Just dorsal versus ventral
    allMeansBtw=( stats.nanmean(dorsBtwCoher), stats.nanmean(ventBtwCoher))
    allSTDBtw=( stats.nanstd(dorsBtwCoher), stats.nanstd(ventBtwCoher))

    # Make bar graph
    title= titleName+ 'by Network for ' +sub+ ' for '+ str(numRuns)+' runs'; labels=( 'Dorsal', 'Ventral')
    makeBarPlots(allMeansWithin, allSTDWithin, allMeansBtw, allSTDBtw, title, labels)

def getNetworkWithin(in_dat, roiIndx):
    m=in_dat.copy()
    #Null the upper triangle, so that you don't get the redundant and
    #the diagonal values:
    idx_null = triu_indices(m.shape[0])
    m[idx_null] = np.nan

    #Extract network values
    withinVals = m[roiIndx,:][:,roiIndx]
    return withinVals


def getNetworkBtw(dataBtw, net1, net2):
    data_b=dataBtw.copy()
    allBtw=data_b[net1,:][:,net2]
    allBtw=stats.nanmean(allBtw)
    return allBtw
