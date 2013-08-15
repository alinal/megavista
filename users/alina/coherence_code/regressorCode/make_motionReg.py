# Preprocessing code
#
# Time series are already motion corrected.
#
#

import subjects_regressors
reload(subjects_regressors) # In case you make changes in there while you analyze
from subjects_regressors import subjects

from rd_afni_make_1D import txt_to_1D

if __name__ == "__main__":

    base_path = '/Volumes/Plata1/DorsalVentral/' # Change this to your path
    fmri_path = base_path + 'fmri/'

    sessionName=['donepazil', 'placebo']
    session=[0,1] # 0= donepazil, 1=placebo
    TR = 2
    allRuns=['fix_nii', 'right_nii', 'left_nii']

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
                    txt_to_1D(txt_file, out_file_base)

