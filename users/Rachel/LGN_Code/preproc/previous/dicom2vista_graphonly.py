#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
#dicom2vista.py

"""This script takes the raw dicoms as they come off the scanner, and
creates the environment that mrVista expects to find, when running
mrInit2

100109 ASR wrote it
""" 

import os
import glob
import numpy as np
import sys
from matplotlib import pyplot as plt

if __name__ == "__main__": 
    #The full path to the session files is a command-line argument: 
    sess_dir = sys.argv[1]
    if sess_dir[-1]=='/': #If a trailing backslash has been input
        sess_dir=sess_dir[:-1]
    sess_name = os.path.split(sess_dir)[1]
    #switch to session directory:
    nifti_dir = sess_dir + '/' + sess_name + '_nifti'
    os.chdir(nifti_dir)
    print os.path.realpath(os.path.curdir)
    dir_list = np.array(os.listdir('.')) 
    motion_list = []
    for file in dir_list:
        if file.endswith('mcf.par'):
            motion_list.append(file)
    dir_list = motion_list
    print dir_list

    #make empty containers to add the motion params to, for plotting purposes:
    #Three rotation params:
    r1 = np.array([]) 
    r2 = np.array([])
    r3 = np.array([])
    #And three translations:
    t1 = np.array([])
    t2 = np.array([])
    t3 = np.array([])

    for this_dir in dir_list: 
        dt = dict(names = ('R1','R2','R3','T1','T2','T3'), 
                  formats = (np.float32,np.float32,np.float32,
                  np.float32,np.float32,np.float32))
        motion_params = np.loadtxt(this_dir,dt)
            
        r1 = np.append(r1,motion_params['R1'])
        r2 = np.append(r2,motion_params['R2'])
        r3 = np.append(r3,motion_params['R3'])
        t1 = np.append(t1,motion_params['T1'])
        t2 = np.append(t2,motion_params['T2'])
        t3 = np.append(t3,motion_params['T3'])

    #After doing all that copy stuff into the right places:
    for this_dir in dir_list:
        print os.path.realpath(os.path.curdir)


    #os.chdir('../../' + sess_name + '_nifti/')
    
    #Plot the motion params: 
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(t1)
    ax1.plot(t2)
    ax1.plot(t3)
    ax1.set_ylabel('Translation (mm)')
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(r1)
    ax2.plot(r2)
    ax2.plot(r3)
    ax2.set_ylabel('Rotation (rad)')
    ax2.set_xlabel('Time (TR)')
    fig.savefig(sess_name + '_motion_params.png')
    os.system
    ('open ' + sess_name + '_motion_params.png')
