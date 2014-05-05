# This code loads each ROI's time series and does the following:
# 1. Linearly detrend using signal.detrend
# 2. Band pass filter the data using boxcar filter in filter analyzer (default filter size)
# 3. Least squares regression of nuisance variables (motion and ROI) using np.linalg.lstsq
# 4. Takes the residuals from the regression
# 5. Averages over all ROI voxels
# 6. Saves out numRun x numROI x TR matrix of data

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
import scipy.linalg

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


if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

    plotAll=0;
    saveROI=1;

    sessionName=['donepazil', 'placebo']
    session=[1] # 0= donepazil, 1=placebo
    TR = 2
    allRuns=['right_nii', 'left_nii', 'fix_nii']
    #allRuns=['fix_nii']
    # TR=2 seconds, 30 TRs in one movie
    segTime=30

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

            # Get ROIs
            #roi_names=np.array(rois)
            roi_names=np.array(rois[subject][sess][1])
            ROI_files=[]
            for roi in rois[subject][sess][1]:
                ROI_files.append(fmri_path+sessName[0]+'/Inplane/ROIs/' +roi +'.mat')

            # Get the coordinates of the ROIs, while accounting for the
            # up-sampling:
            ROI_coords = [tsv.upsample_coords(tsv.getROIcoords(f),up_samp) for f in ROI_files]


            nifti_path = fmri_path +sessName[0] + '/%s_nifti/' % sessName[0]
            reg_path=fmri_path+sessName[0]+'/regressors/'

            # Add filtering
            filterType='boxcar'

            #Go through each run and save out ROIs as nifti file
            for runName in allRuns:
                print 'Analyzing ' + runName
                # Initialize lists for each condition:
                t_fix = []
                saveFile=base_path+ 'fmri/Results/timeseries/'+subject+sessionName[sess]+'_'+runName+'_%sROIts_%sReg_stc.pck' % (len(roi_names), len(nuisReg))
                for this_run in sessName[1][runName]:
                    run_rois=[]
                    # Load stc nifti
                    allData=load_nii(nifti_path+this_run[:-7]+'_stc.nii.gz', ROI_coords, TR, normalize='percent', average=False, verbose=True)
                    regMatrix=[]

                    # Load the nuisance variables into a matrix
                    for reg in nuisReg:
                        regFile=reg_path+this_run[:5]+'_'+reg
                        regMatrix.append(np.loadtxt(regFile))
                        print 'Regressor ' + reg
                    # Convert to array
                    regArray=np.array(regMatrix).transpose()

                    # Go through each ROI
                    for jj in range(len(ROI_coords)):
                        print 'Analyzing '+ roi_names[jj]
                        roiData=[]; ts_roidt=[]; ts_Box=[]; ts_AvgBox=[];
                        # Get time series for each ROI
                        roiData=allData[jj]

                        # Linearly detrend within each voxel
                        ts_roidt=signal.detrend(roiData.data, axis=1)

                        # Band pass filter the data using boxcar filter
                        ts_Box=bp_data(ts_roidt, TR, f_ub, f_lb)

                        # Average over all voxels
                        ts_AvgBox=np.mean(ts_Box.data, 0)

                        if plotAll:
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

                        # Regression steps
                        tsArray=np.array(ts_Box.data)
                        residMatrix=[]

                        # Do multiple regression using least squares.
                        # Regress within each voxel
                        for vox in range(tsArray.shape[0]):
                            b_weight=[]
                            b_weight=np.linalg.lstsq(regArray,tsArray[vox])[0]
                            residMatrix.append(tsArray[vox]-np.dot(regArray,b_weight))

                        # Make residual array
                        residArray=np.array(residMatrix)

                        # Get ROI average
                        roiAvg=[]; roiAvg=np.mean(residArray, 0)
                        run_rois.append(ts.TimeSeries(roiAvg, sampling_interval=TR))
                    t_fix.append(run_rois)

                if saveROI:
                    # Save time series
                    file=open(saveFile, 'w') # write mode
                    # First file loaded is TS files
                    pickle.dump(t_fix, file)
                    # Second file loaded is ROI names
                    pickle.dump(roi_names, file)
                    # Subject and type
                    pickle.dump(subject+sessionName[sess], file)
                    # Third file is run names
                    pickle.dump(sessName[1][runName], file)
                    # Session and run
                    pickle.dump(sessName[0]+runName, file)
                    file.close()
                    print 'Saving subject dictionaries in %s.' % saveFile
