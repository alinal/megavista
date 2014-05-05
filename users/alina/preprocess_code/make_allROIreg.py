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

import subjects_roi_regressors
reload(subjects_roi_regressors) # In case you make changes in there while you analyze
from subjects_roi_regressors import subjects, rois
from preproc_filter import bp_data

if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

    sessionName=['donepazil', 'placebo']
    session=[1] # 0= donepazil, 1=placebo
    TR = 2
    allRuns=['right_nii', 'left_nii', 'fix_nii']
    allRuns=['fix_nii']

    # The pass band is f_lb <-> f_ub.
    # Also, see: http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
    f_ub = 0.15
    f_lb = 0.01

    # The upsample factor between the Inplane and the Gray:
    # Inplane Voxels: .867 x .867 x 3.3, Functional voxels: 3 x 3 x 3.3
    up_samp = [3.4595,3.4595,1.0000]


    for subject in subjects:
        # Get session
        for sess in session:

            sessName = subjects[subject][sess]
            print 'Analyzing ' + sessName[0]

            # Get ROIs
            roi_names=np.array(rois[subject][sess][1])
            ROI_files=[]
            for roi in rois[subject][sess][1]:
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

            #Go through each run and save out ROIs as nifti file
            # Band pass filtering and detrending are done on ROI average timeseries
            for runName in allRuns:
                print 'Analyzing ' + runName
                for this_fix in sessName[1][runName]:
                    t_all=[]
                    # Load the time series and average over ROI
                    t_all=load_nii(nifti_path+this_fix[:-4]+'_stc.nii.gz', ROI_coords,TR, normalize='percent', average=True, verbose=True)

                    for roiNum in range(len(roi_names)):
                        print 'Analyzing '+ roi_names[roiNum]
                        ts_roi=[]; ts_roidt=[]; ts_Box=[];
                        ts_roidv=[]; ts_roidv_dt=[]; ts_roidv_dtBox=[];

                        # Get each time series (1 x TRs)
                        ts_roi=t_all[roiNum].data
                        # Linearly detrend each ROI
                        ts_roidt=signal.detrend(ts_roi, axis=0)

                        # Band pass filter the data using boxcar filter
                        ts_Box=bp_data(ts_roidt, TR, f_ub, f_lb)
                        1/0

                        # Get the derivative
                        ts_roidv=np.diff(ts_roi)
                        # Make it the same size as the original
                        ts_roidv=np.insert(ts_roidv, 0, 0)
                        # Detrend the derivative
                        ts_roidv_dt=signal.detrend(ts_roidv)
                        #Bandpass the derivative
                        ts_roidv_dtBox=bp_data(ts_roidv_dt, TR, f_ub, f_lb)

                        # Plot TS results
                        #plt.figure(); plt.plot(ts_roi); plt.plot(ts_roidt); plt.plot(ts_Box) ;
                        #plt.legend(('Original TS', 'Linearly Filtered TS', 'Bandpass filtered'))

                        # Plot frequencies
                        S_boxcar=SpectralAnalyzer(ts.TimeSeries(ts_Box, sampling_interval=TR))
                        S_original=SpectralAnalyzer(ts.TimeSeries(t_all[roiNum], sampling_interval=TR))
                        S_dt=SpectralAnalyzer(ts.TimeSeries(ts_roidt, sampling_interval=TR))

                        #fig03 = plt.figure()
                        #ax03 = fig03.add_subplot(1, 1, 1)
                        #ax03.plot(S_original.spectrum_multi_taper[0], S_original.spectrum_multi_taper[1], label='Original')
                        #ax03.plot(S_dt.spectrum_multi_taper[0], S_dt.spectrum_multi_taper[1], label='Detrended')
                        #ax03.plot(S_boxcar.spectrum_multi_taper[0], S_boxcar.spectrum_multi_taper[1], label='Boxcar')
                        #ax03.legend()
                        1/0
                        # Save nuisance time series
                        out_file=save_path+this_fix[:-8]+'_'+roi_names[roiNum]+'_stc_avgFt.1D'
                        np.savetxt(out_file, ts_Box)
                        out_file=save_path+this_fix[:-8]+'_'+roi_names[roiNum]+'_stc_avgFt_deriv.1D'
                        np.savetxt(out_file, ts_roidv_dtBox)


