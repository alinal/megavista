#!/usr/bin/python

import os, string
import numpy as np

## DEFINITIONS
def gzip_files(topDir, afni_file,localizerfile):
    """ Script which gzips files not currently being used """

    print 'GZIPPING FILES...'

    command1 = "gzip %s*BRIK" % (topDir)
    command2 = "gzip %s*HEAD" % (topDir)
    command3 = "gzip %s%s/pb00*" % (topDir,afni_file)
    command4 = "gzip %s%s/errts*" % (topDir,afni_file)
    command5 = "gzip %s%s/fitts*" % (topDir,afni_file)
    command6 = "gzip %s%s/*al*" % (topDir,afni_file)
    command7 = "gzip %s%s/*ns*" % (topDir,afni_file)
    command8 = "gzip %s%s/all_runs*" % (topDir,afni_file)
    command9 = "gzip %s%s/*mask*" % (topDir,afni_file)
    command10 = "gzip %s*nii" % (topDir)
    command11 = "gzip %s%s/pb01*" % (topDir,afni_file)
    command12 = "gzip %s%s/pb03*" % (topDir,afni_file)
    command13 = "gzip %s%s/pb02*" % (topDir,afni_file)

    os.system(command1)
    os.system(command2)
    os.system(command3)
    os.system(command4)
    os.system(command5)
    os.system(command6)
    os.system(command7)
    os.system(command8)
    #os.system(command9)
    os.system(command10)
    os.system(command11)
    try:
        os.system(command12)
    except:
        print 'Could not gzip pb03 files'

    if localizerfile == 'facescene_loc':
        os.system(command13)
        

##INPUTS
#subs = ['S10D1','S10D2','S11D1','S11D2']
#['SS2','SS3']
#['S07D1','S07D2','S08D1','S08D2','S10D1','S10D2','S11D1','S11D2']
#['EP1','RM1','RM2','SS2','SS3','TC4','TC5'] #[EN3] #['EP2'] ['AM1','AM3','EN2']
#subs = ['AM1','AM3','EN2','EN3','EP1','EP2']
#subs = ['RM1','RM2','SS2','SS3','TC4','TC5']
subs = ['S07D1','S07D2','S08D1','S08D2','S10D1','S10D2','S11D1','S11D2']

#What to do
epi_space = False #always false now...
make_script = False #can't be true WITH run_script
run_script = True  #***need to change duration of basis functions on script before running
gzip_all = False

##MAIN SCRIPT
project = 'FaceSpace_loc'
filelist = ['faceattend']

for localizerfile in filelist:
    for sub in subs:

        #code for scan type
        if sub == 'AM1' or sub == 'AM3' or sub == 'EN2':
            scan_type = 'highSNR'
        else:
            scan_type = 'highRES'


        #quick sanity check:
        if run_script and make_script:
            print 'NEVER make script while running it!'
            1/0


        #filename
        topDir = '/home/despo/cgratton/data/FaceSpace_loc/Data/'+sub+'/' +localizerfile+ '_' + scan_type + '/Analysis/'

        if epi_space is False:
            if localizerfile == 'facescene_loc':
                #command = 'bash do_afni_proc_localizer.txt %s %s %s' % (sub,localizerfile,scan_type)
                1/0
                
            elif localizerfile == 'facespace_loc':
                #command = 'bash do_afni_proc_facespace_respprd.txt %s %s %s' % (sub,localizerfile,scan_type)
                1/0
                
            elif localizerfile == 'faceattend':
                #command = 'bash do_afni_proc_faceattend_respprd.txt %s %s %s' % (sub,localizerfile,scan_type)

                if sub == 'TC4' or sub == 'TC5':
                    command = 'bash do_afni_proc_faceattend_STC_blur.txt %s %s %s' % (sub,localizerfile,scan_type)
                else:
                    command = 'bash do_afni_proc_faceattend_respprd_STC_blur.txt %s %s %s' % (sub,localizerfile,scan_type)

            if sub == 'TC4' or sub == 'TC5':
                afni_file = 'afni_proc_STC_blur'
            else:
                afni_file = 'afni_proc_respprd_STC_blur'




        if make_script:
            os.system(command)




        if run_script:
            command = """tcsh -xef /home/despo/cgratton/data/FaceSpace_loc/Scripts/%s_%s_%s.%s.script""" % (sub,localizerfile,scan_type,afni_file)

            print command
            os.system(command)



        if gzip_all:
            gzip_files(topDir,afni_file,localizerfile)
