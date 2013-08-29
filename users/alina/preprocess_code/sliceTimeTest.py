# Changed on 1/28/11
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


import subjects
reload(subjects) # In case you make changes in there while you analyze
from subjects import subjects, rois

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


if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'
    sessionName=['donepazil', 'placebo']
    session=1 # 0= donepazil, 1=placebo
    TR = 2
    allRuns=['fix_nii']

    # The upsample factor between the Inplane and the Gray:
    # Inplane Voxels: .867 x .867 x 3.3, Functional voxels: 3 x 3 x 3.3
    up_samp = [3.4595,3.4595,1.0000]

    for subject in subjects:
            # len(subjects[subject])= number of session per subject
            # len(subjects[subject][0][1])= number of different types of runs
            # len(subjects[subject][1][1]['fix_nii'])= number of nifti files for that session

            # Get session
            sess = subjects[subject][session]

            # Get ROIs
            roi_names=np.array(rois)
            ROI_files=[]
            for roi in rois:
                ROI_files.append(fmri_path+sess[0]+'/Inplane/ROIs/' +roi +'.mat')

            # Get the coordinates of the ROIs, while accounting for the
            # up-sampling:
            ROI_coords = [tsv.upsample_coords(tsv.getROIcoords(f),up_samp) for f in ROI_files]

             # Initialize lists for each behavioral condition:
            t_fix = []
            nifti_path =fmri_path+sess[0]+'/%s_nifti/' % sess[0]


            niftiOrig=load_nii(nifti_path+'epi04_mcf.nii.gz', ROI_coords, TR, average='True')
            niftiSTc=load_nii(nifti_path+'epi04_mcf_stc.nii.gz', ROI_coords, TR, average='True')

            roiNum=13;
            roi1Orig=niftiOrig[roiNum]
            roi1ST=niftiSTc[roiNum]

            plt.figure(); plt.plot(roi1Orig.data); plt.plot(roi1ST.data); plt.legend(('Original', 'Slice Time'));
            plt.title(rois[roiNum])

            1/0
