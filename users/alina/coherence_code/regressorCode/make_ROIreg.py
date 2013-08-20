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

from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer, SpectralAnalyzer, FilterAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr

from scipy import signal

import subjects_regressors
reload(subjects_regressors) # In case you make changes in there while you analyze
from subjects_regressors import subjects, rois

if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

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
            nifti_path = fmri_path +sessName[0] + '/%s_nifti/' % sessName[0]
            save_path=fmri_path+sessName[0]+ '/regressors/'

            # Add filtering
            filterType='boxcar'

            #Go through each run and save out motion parameters as 1D text file
            for runName in allRuns:
                for this_fix in sessName[1][runName]:
                    t_all=[]
                    t_all=load_nii(nifti_path+this_fix[:-4]+'.nii.gz', ROI_coords,TR, normalize='percent', average=False, verbose=True)
                    for roiNum in range(len(rois)):
                        ts_roi=[]
                        ts_roidt=[]
                        ts_roidtAvg=[]
                        ts_roidtAvgConv=[]

                        # Get each time series (voxels x TRs)
                        ts_roi=t_all[roiNum].data

                        # Linearly detrend within each voxel
                        ts_roidt=signal.detrend(ts_roi, axis=1)

                        # Average across all voxels within ROI
                        ts_roidtAvg=np.mean(ts_roidt, 0)

                        # Band pass filter the data
                        ts_roidtAvgConv=ts.TimeSeries(ts_roidtAvg, sampling_interval=TR)
                        ts_roidtAvgAnal=FilterAnalyzer(ts_roidtAvgConv, lb=f_lb, ub=f_ub)
                        ts_AvgBox=ts_roidtAvgAnal.filtered_boxcar

                        # Plot TS results
                        1/0
                        origTS=np.mean(ts_roi, 0)
                        plt.figure(); plt.plot(origTS); plt.plot(ts_roidtAvg); plt.plot(ts_AvgBox.data) ;
                        plt.title('Blue= Original TS; Green= Linearly Filtered TS; Red= bandpass filtered')

                        # Plot frequencies
                        S_boxcar=SpectralAnalyzer(ts_AvgBox)
                        S_original=SpectralAnalyzer(ts.TimeSeries(np.mean(t_all[roiNum], 1), sampling_interval=TR))
                        S_dt=SpectralAnalyzer(ts_roidtAvgConv)

                        fig03 = plt.figure()
                        ax03 = fig03.add_subplot(1, 1, 1)

                        ax03.plot(S_original.spectrum_multi_taper[0], S_original.spectrum_multi_taper[1], label='Original')

                        ax03.plot(S_dt.spectrum_multi_taper[0], S_dt.spectrum_multi_taper[1], label='Detrended')

                        ax03.plot(S_boxcar.spectrum_multi_taper[0], S_boxcar.spectrum_multi_taper[1], label='Boxcar')

                        ax03.legend()

                        # Save nuisance time series
                        out_file=save_path+rois[roiNum]+'_'+this_fix[:-8]+'.1D'
                        np.savetxt(out_file, ts_AvgBox.data)


