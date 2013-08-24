
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pickle
import datetime
import copy

import vista_utils as tsv # Get it at: https://github.com/arokem/vista_utils
from nitime.fmri.io import time_series_from_file as load_nii
import nitime.timeseries as ts
import nitime.viz as viz

from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer, SpectralAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

def reshapeTS(t_fix):
    # TR=2 seconds, 30 TRs in one movie
    segTime=30
    # Change to an array (numSess, numROIs, numTime points)
    t_fixArray=np.array(t_fix)
    t_fixArrayTP=np.transpose(t_fixArray, (1,0,2))
    shapeTS=t_fixArrayTP.shape
    numRuns=shapeTS[2]/segTime
    # This returns rois x runs x TS with runs collapsed by segTime
    allROIS=np.reshape(t_fixArrayTP, [shapeTS[0], shapeTS[1]*numRuns, segTime])
    return allROIS


if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    subFiles=['CGplacebo_right_nii_43ROIts.pck', 'CGplacebo_left_nii_43ROIts.pck' ]
    normalizeByMean=1
    plot=0

    # The pass band is f_lb <-> f_ub.
    # Also, see: http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
    f_ub = 0.15
    f_lb = 0.01

    #It depends on your frequency bin resolution needs. delta_freq = sampling_rate/NFFT
    #So, say your sampleing rate is 1024 samples/sec and NFFT is 256. Then delta_freq = 4 Hz.
    NFFT=16 # 32 for 60 TRs, 1/64= freq limit lower, .25 hertz is upper limit (1/2 of sampling rate) Nyquist freq
    n_overlap=8
    TR = 2

    filterType='boxcar'

    for subject in range(len(subFiles)):
        loadFile=base_path+ 'fmri/Results/timeseries/'+subFiles[subject]
        date=str(datetime.date.today())

        # set up dictionaries to store results
        corr_all=dict()
        coh_all = dict()
        hierarch=dict()

        # load time series
        print 'Loading %s.' % loadFile
        file=open(loadFile, 'r')
        # First file is TS
        t_fix=pickle.load(file)
        # Next file is roi names
        roiNames=pickle.load(file)
        # Subject and type
        subName=pickle.load(file)
        runNums=pickle.load(file)
        subInfo=pickle.load(file)
        file.close()

        # reshape ROI matrix
        allROIS=reshapeTS(t_fix)
        numRuns=allROIS.shape[1]

        corr_all[subName] = np.zeros((numRuns,len(roiNames),len(roiNames))) * np.nan
        coh_all[subName] = np.zeros((numRuns,len(roiNames),len(roiNames))) * np.nan
        hierarch[subName]=np.zeros((numRuns, len(roiNames)*(len(roiNames)-1)/2)) * np.nan

        # Average over all the runs, get an ROI by TS array (TS==averages), TS= 30 TRs long (TR=2 S)
        allROISorig=copy.deepcopy(allROIS)
        AvgRuns=np.mean(allROIS, axis=1)

        # Subtract out average
        if normalizeByMean:
            for numTS in range(numRuns):
                allROIS[:,numTS,:]=allROIS[:,numTS,:]-AvgRuns
                # plt.figure()
                #Examples time series subtraction, one ROI
                #plt.plot(allROIS[0,numTS,:], color='r', label='newTS')
                #plt.plot(allROISorig[0,numTS,:], color='b', label='oldTS')
                #plt.plot(AvgRuns[0,:], color='g', label='avgTS')
                #plt.legend(loc='best'); plt.title(roi_names[0]);
                #plt.show()

      # Get roi correlations and coherence
        for runNum in range(allROIS.shape[1]):
            #need to load timeseries by run
            fixTS=ts.TimeSeries(allROIS[:,runNum,:], sampling_interval=TR)
            fixTS.metadata['roi'] = roiNames

            # Get plot and correlations
            C=CorrelationAnalyzer(fixTS)
            #fig01 = drawmatrix_channels(C.corrcoef, roi_names, size=[10., 10.], color_anchor=0,  title='Correlation Results Run %i' % runNum)
            #plt.show()

            # Save correlation
            corr_all[subName][runNum]=C.corrcoef

            # Get coherence
            Coh = CoherenceAnalyzer(fixTS)

            Coh.method['NFFT'] = NFFT
            Coh.method['n_overlap']=n_overlap

            # Get the index for the frequencies inside the ub and lb
            freq_idx = np.where((Coh.frequencies > f_lb) * (Coh.frequencies < f_ub))[0]

            # Extract coherence
            # Coher[0]= correlations for first ROI in list with others
            coher = np.mean(Coh.coherence[:, :, freq_idx], -1)  # Averaging on the last dimension
            #fig03 = drawmatrix_channels(coher, roi_names, size=[10., 10.], color_anchor=0, title='Coherence Results Run %i' % runNum)
            # Save coherence (coher is the average of the coherence over the specified frequency)
            coh_all[subName][runNum]=coher

            # Get hierarchical clustering distances
            hierarch[subName][runNum]=pdist(fixTS.data)

            if plot:
                #For debugging, lets look at some of the spectra
                S_original = SpectralAnalyzer(fixTS)
                plt.figure()
                plt.plot(S_original.psd[0],S_original.psd[1][0],label='Welch PSD')
                plt.plot(S_original.spectrum_fourier[0],S_original.spectrum_fourier[1][0],label='FFT')
                plt.plot(S_original.periodogram[0],S_original.periodogram[1][0],label='Periodogram')
                plt.plot(S_original.spectrum_multi_taper[0],S_original.spectrum_multi_taper[1][0],label='Multi-taper')
                plt.xlabel('Frequency (Hz)')
                plt.ylabel('Power')
                plt.ylim([-5, 5])
                plt.legend()
                plt.title('Run number '+str(runNum+1)+', '+filterType)
                plt.show()

        saveFile=base_path+'fmri/Results/correlation/'+subFiles[subject][:-10]+'_corrVals_wGM_hierarch.pck'

        file=open(saveFile, 'w') # write mode
        # First file loaded is coherence
        pickle.dump(coh_all, file)
        # Second file loaded is correlation
        pickle.dump(corr_all, file)
        # Save roi names
        pickle.dump(roiNames, file)
        # save subjects
        pickle.dump(subFiles, file)
        # Save hierarchical clustering distances
        pickle.dump(hierarch, file)
        file.close()
        print 'Saving coherence and correlation dictionaries in %s.' % saveFile





