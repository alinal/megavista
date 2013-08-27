# Preprocessing code
#
# Time series are already motion corrected.
#
#
# Combined with RD_txt_to_1D code to save txt
#   """ txt_file is a text file with array-like data
#    the columns of txt_file will be saved as separate .1D files
#    out_file_base is the prefix for the saved files
#    the column number (1...ncolumns) will be appended to the prefix
#    """

import sys
import os

import matplotlib.pyplot as plt

import subjects_regressors
reload(subjects_regressors) # In case you make changes in there while you analyze
from subjects_regressors import subjects
import numpy as np

import nitime.timeseries as ts
from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer, SpectralAnalyzer, FilterAnalyzer

from preproc_filter import bp_data

if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

    allRuns=['fix_nii', 'right_nii', 'left_nii']
    sessionName=['donepazil', 'placebo']
    session=[0,1] # 0= donepazil, 1=placebo
    TR = 2

    # The pass band is f_lb <-> f_ub.
    # Also, see: http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
    f_ub = 0.15
    f_lb = 0.01

    for subject in subjects:

        # Get session
        for sess in session:
            sessName = subjects[subject][sess]

            # Initialize lists for each behavioral condition:
            t_fix = []
            nifti_path = fmri_path +sessName[0] + '/%s_nifti/' % sessName[0]
            save_path=fmri_path+sessName[0]+ '/regressors/'

            #Go through each run and save out motion parameters as 1D text file
            for runName in allRuns:
                for this_fix in sessName[1][runName]:
                    txt_file=nifti_path+this_fix
                    out_file_base=save_path+this_fix[:-4]

                    # read in par file
                    data = np.loadtxt(txt_file)

                    # save the columns of data as separate .1D files
                    for i in xrange(data.shape[1]):
                        mt_deriv=[];

                        # Get derivative
                        mt_deriv=np.diff(data[:,i])
                        # Make it the same size as the original
                        mt_deriv=np.insert(mt_deriv, 0, 0)

                        # Run bandpass filter
                        mt_bp=bp_data(data[:,i], TR, f_ub, f_lb)
                        mt_bp_deriv=bp_data(mt_deriv, TR, f_ub, f_lb)

                        # Save file
                        out_file = '{0}{1}_bp.1D'.format(out_file_base, i+1)
                        np.savetxt(out_file, mt_bp)

                        # Save file
                        out_file = '{0}{1}_bp_deriv.1D'.format(out_file_base, i+1)
                        np.savetxt(out_file, mt_bp_deriv)


