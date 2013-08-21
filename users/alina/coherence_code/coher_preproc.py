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

from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer, SpectralAnalyzer,  FilterAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr
from preproc_filter import bp_data
from scipy import signal

import subjects
reload(subjects) # In case you make changes in there while you analyze
from subjects import subjects, rois, nuisReg

def display_vox(tseries,vox_idx,fig=None):
    """
    Display the voxel time-series
    """
    if fig is None:
        fig = plt.figure()

    vox_tseries = ts.TimeSeries(tseries.data[vox_idx],sampling_interval=TR)

    fig = viz.plot_tseries(vox_tseries,fig)
    fig = viz.plot_tseries(ts.TimeSeries(np.mean(vox_tseries.data,0),
                                         sampling_interval=TR),
                           yerror=ts.TimeSeries(stats.sem(vox_tseries.data,0),
                                                sampling_interval=TR),fig=fig,
                           error_alpha = 0.3,ylabel='% signal change',
                           linewidth=4,
                           color='r')
    return fig

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

def loadNuisance(files):


if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

    plotAll=1;

    sessionName=['donepazil', 'placebo']
    session=[0,1] # 0= donepazil, 1=placebo
    TR = 2
    #allRuns=['fix_nii', 'right_nii', 'left_nii']
    allRuns=['fix_nii']

    # The pass band is f_lb <-> f_ub.
    # Also, see: http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
    f_ub = 0.15
    f_lb = 0.01

    #It depends on your frequency bin resolution needs. delta_freq = sampling_rate/NFFT
    #So, say your sampleing rate is 1024 samples/sec and NFFT is 256. Then delta_freq = 4 Hz.
    NFFT=16 # 32 for 60 TRs, 1/64= freq limit lower, .25 hertz is upper limit (1/2 of sampling rate) Nyquist freq
    n_overlap=8

    # The upsample factor between the Inplane and the Gray:
    # Inplane Voxels: .867 x .867 x 3.3, Functional voxels: 3 x 3 x 3.3
    up_samp = [3.4595,3.4595,1.0000]


    for subject in subjects:
        # Get session
        for sess in session:
            sessName = subjects[subject][sess]

            # Get ROIs
            roi_names=np.array(rois)
            ROI_files=[]
            for roi in rois:
                ROI_files.append(fmri_path+sessName[0]+'/Inplane/ROIs/' +roi +'.mat')

            # Get the coordinates of the ROIs, while accounting for the
            # up-sampling:
            ROI_coords = [tsv.upsample_coords(tsv.getROIcoords(f),up_samp)
                           for f in ROI_files]

            # Initialize lists for each behavioral condition:
            t_fix = []
            nifti_path = fmri_path +sessName[0] + '/%s_nifti/' % sessName[0]
            reg_path=fmri_path+sessName[0]+'/regressors/'

            # Add filtering
            filterType='boxcar'

            #Go through each run and save out ROIs as nifti file
            for runName in allRuns:
                for this_fix in sessName[1][runName]:
                    for jj in range(len(ROI_coords)):
                        roiData=[]; ts_roidt=[]; ts_Box=[]; ts_AvgBox=[];
                        1/0
                        # Load time series for each ROI
                        roiData=load_nii(nifti_path+this_fix, ROI_coords[jj], TR, normalize='percent', average=False, verbose=True)

                        # Linearly detrend within each voxel
                        ts_roidt=signal.detrend(roiData.data, axis=1)

                        # Band pass filter the data using boxcar filter
                        ts_Box=bp_data(ts_roidt, TR, f_ub, f_lb)

                        # Average over all voxels
                        ts_AvgBox=np.mean(ts_Box.data, 0)

                        if plotAll
                            # Plot time series
                            origTS=np.mean(roiData.data, 0)
                            plt.figure(); plt.plot(origTS);  plt.plot(ts_AvgBox) ;
                            plt.legend(('Original TS', 'Bandpass filtered'))

                            # Plot frequencies
                            S_boxcar=SpectralAnalyzer(ts.TimeSeries(ts_AvgBox, sampling_interval=TR))
                            S_original=SpectralAnalyzer(ts.TimeSeries(np.mean(roiData.data, 0), sampling_interval=TR))

                            fig03 = plt.figure()
                            ax03 = fig03.add_subplot(1, 1, 1)

                            ax03.plot(S_original.spectrum_multi_taper[0], S_original.spectrum_multi_taper[1], label='Original')

                            ax03.plot(S_boxcar.spectrum_multi_taper[0], S_boxcar.spectrum_multi_taper[1], label='Boxcar')

                            ax03.legend()

                        # Do regression
                        regMatrix=[]
                        for reg in nuisReg
                            regMatrix.append(np.loadtxt(nifti_))

                        # Load the nuisance variables into a matrix

                        1/0
                    t_fix.append(load_nii(nifti_path+this_fix, ROI_coords, TR, normalize='percent', average=False, verbose=True))

            #for roiNum in range(len(rois)):
            # Get each time series (voxels x TRs)