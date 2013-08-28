'''
 Data structure used in all analysis programs, containing the lists of scan
 numbers in each condition for each subject
 mcf= already motion corrected

 Subjects are watching a movie: fix= task at the center, right= attention right, contrast decrement


 For each subject, the first session is the donepezil session
 Each run has 6 movie clips x 60 seconds each
'''


# Leave out: 'r_IOG_p3' 'r_pSTS_p3', r_LOf_p3,  'R_MT_al_.5_0.25',
#DCA: no MT
# CHT: no rIPS5
# 'R_V1_0.25', 'R_V1_0.25',
# WC: no 'r_pFus_p3', 'r_mFus_p3', 'r_PPA_p4',
# Get invalid index error with WC
# SS: no 'R_IPS3_0.25', 'R_IPS4_0.25', 'R_IPS5_0.25'
# CG: no l_pFus

rois=['rVent_4mm', 'lVent_4mm', 'rWM_6mm', 'lWM_6mm', 'gray_al'];

subjects = { 'DCA':[['DCA042511',dict(loc_nii =['epi01_mcf.par',
                                       'epi11_mcf.par'],
                             fix_nii = ['epi04_mcf.par',
                                        'epi07_mcf.par',
                                        'epi10_mcf.par'],
                             left_nii =  ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'],
                             right_nii = ['epi02_mcf.par',
                                          'epi05_mcf.par',
                                          'epi08_mcf.par'])],
                ['DCA041111', dict(loc_nii =['epi01_mcf.par',
                                             'epi11_mcf.par'],
                                   fix_nii = ['epi04_mcf.par',
                                              'epi07_mcf.par',
                                              'epi10_mcf.par'],
                                   left_nii = ['epi03_mcf.par',
                                                'epi06_mcf.par',
                                                'epi09_mcf.par'],
                                   right_nii = ['epi02_mcf.par',
                                                'epi05_mcf.par',
                                                'epi08_mcf.par'])]]}


'''

'CHT':[['CHT042111', dict(loc_nii =['epi01_mcf.par',
                                       'epi11_mcf.par'],
                             right_nii = ['epi04_mcf.par',
                                        'epi07_mcf.par',
                                        'epi10_mcf.par'],
                             left_nii =  ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'],
                             fix_nii = ['epi02_mcf.par',
                                        'epi05_mcf.par',
                                        'epi08_mcf.par'])],
            ['CHT051911',dict(loc_nii =['epi01_mcf.par',
                                        'epi11_mcf.par'],
                             right_nii = ['epi04_mcf.par',
                                          'epi07_mcf.par',
                                          'epi10_mcf.par'],
                             left_nii = ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'],
                             fix_nii = ['epi02_mcf.par',
                                        'epi05_mcf.par',
                                        'epi08_mcf.par'])]]}

'CG':[['CG011611', dict(loc_nii =['epi01_mcf.par',
                                      'epi10_mcf.par'],
                            fix_nii = ['epi04_mcf.par',
                                       'epi07_mcf.par'],
                            right_nii =  ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'],
                            left_nii = ['epi02_mcf.par',
                                        'epi05_mcf.par',
                                        'epi08_mcf.par'])],
          ['CG020611', dict(loc_nii = ['epi01_mcf.par',
                                       'epi11_mcf.par'],
                            fix_nii = ['epi04_mcf.par',
                                       'epi07_mcf.par',
                                       'epi10_mcf.par'],
                            right_nii = ['epi03_mcf.par',
                                         'epi06_mcf.par',
                                         'epi09_mcf.par'],
                            left_nii =  ['epi02_mcf.par',
                                         'epi05_mcf.par',
                                         'epi08_mcf.par'])]]}




    'WC':[['WC031911', dict(loc_nii =['epi01_mcf.par',
                                       'epi11_mcf.par'],
                             left_nii = ['epi02_mcf.par',
                                        'epi05_mcf.par',
                                        'epi08_mcf.par'],
                             fix_nii =  ['epi04_mcf.par',
                                          'epi07_mcf.par',
                                          'epi10_mcf.par'],
                             right_nii = ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'])],
           ['WC030311', dict(loc_nii =['epi01_mcf.par',
                                       'epi11_mcf.par'],
                             left_nii = ['epi02_mcf.par',
                                        'epi05_mcf.par',
                                        'epi08_mcf.par'],
                             fix_nii =  ['epi04_mcf.par',
                                          'epi07_mcf.par',
                                          'epi10_mcf.par'],
                             right_nii = ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'])]],


    'SS':[['SS012311', dict(loc_nii =['epi01_mcf.par',
                                       'epi11_mcf.par'],
                             fix_nii = ['epi04_mcf.par',
                                        'epi07_mcf.par',
                                        'epi10_mcf.par'],
                             left_nii =  ['epi03_mcf.par',
                                          'epi06_mcf.par',
                                          'epi09_mcf.par'],
                             right_nii = ['epi02_mcf.par',
                                          'epi05_mcf.par',
                                          'epi08_mcf.par'])],
           ['SS011011',dict(loc_nii = ['epi01_mcf.par',
                                       'epi10_mcf.par'],
                            fix_nii = ['epi04_mcf.par',
                                       'epi07_mcf.par'],
                            left_nii = ['epi03_mcf.par',
                                        'epi06_mcf.par',
                                        'epi09_mcf.par'],
                            right_nii =  ['epi02_mcf.par',
                                          'epi05_mcf.par',
                                          'epi08_mcf.par'])]]
                                        }

    # 'MO':[['MO040411', dict(loc_nii =['epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'],
    #                         fix_nii = ['epi02_mcf.nii.gz',
    #                                    'epi05_mcf.nii.gz',
    #                                    'epi08_mcf.nii.gz'],
    #                         left_nii =  ['epi04_mcf.nii.gz',
    #                                      'epi07_mcf.nii.gz',
    #                                      'epi10_mcf.nii.gz'],
    #                         right_nii = ['epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'])],
    #       ['MO041811', dict(loc_nii =['epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'],
    #                         fix_nii = ['epi02_mcf.nii.gz',
    #                                    'epi05_mcf.nii.gz',
    #                                    'epi08_mcf.nii.gz'],
    #                         left_nii =  ['epi04_mcf.nii.gz',
    #                                      'epi07_mcf.nii.gz',
    #                                      'epi10_mcf.nii.gz'],
    #                         right_nii = ['epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'])]],

             # 'DCA':[['DCA042511',dict(loc_nii =['epi01_mcf.nii.gz',
             #                           'epi11_mcf.nii.gz'],
             #                 fix_nii = ['epi04_mcf.nii.gz',
             #                            'epi07_mcf.nii.gz',
             #                            'epi10_mcf.nii.gz'],
             #                 left_nii =  ['epi03_mcf.nii.gz',
             #                              'epi06_mcf.nii.gz',
             #                              'epi09_mcf.nii.gz'],
             #                 right_nii = ['epi02_mcf.nii.gz',
             #                              'epi05_mcf.nii.gz',
             #                              'epi08_mcf.nii.gz'])],
             #    ['DCA041111', dict(loc_nii =['epi01_mcf.nii.gz',
             #                                 'epi11_mcf.nii.gz'],
             #                       fix_nii = ['epi04_mcf.nii.gz',
             #                                  'epi07_mcf.nii.gz',
             #                                  'epi10_mcf.nii.gz'],
             #                       left_nii = ['epi03_mcf.nii.gz',
             #                                    'epi06_mcf.nii.gz',
             #                                    'epi09_mcf.nii.gz'],
             #                       right_nii = ['epi02_mcf.nii.gz',
             #                                    'epi05_mcf.nii.gz',
             #                                    'epi08_mcf.nii.gz'])]],



    # 'JM':[['JM021811', dict(loc_nii =['epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'],
    #                         left_nii = ['epi04_mcf.nii.gz',
    #                                    'epi07_mcf.nii.gz',
    #                                    'epi10_mcf.nii.gz'],
    #                         right_nii =  ['epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'],
    #                         fix_nii = ['epi02_mcf.nii.gz',
    #                                      'epi05_mcf.nii.gz',
    #                                      'epi08_mcf.nii.gz'])],
    #       ['JM030711',dict(loc_nii =['epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'],
    #                         left_nii = ['epi04_mcf.nii.gz',
    #                                    'epi07_mcf.nii.gz',
    #                                    'epi10_mcf.nii.gz'],
    #                         right_nii =  ['epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'],
    #                         fix_nii = ['epi02_mcf.nii.gz',
    #                                      'epi05_mcf.nii.gz',
    #                                      'epi08_mcf.nii.gz'])]],

    # 'MAS':[['MAS030311', dict(loc_nii =['epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'],
    #                         fix_nii = ['epi04_mcf.nii.gz',
    #                                    'epi07_mcf.nii.gz',
    #                                    'epi10_mcf.nii.gz'],
    #                         left_nii =  ['epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'],
    #                         right_nii = ['epi02_mcf.nii.gz',
    #                                      'epi05_mcf.nii.gz',
    #                                      'epi08_mcf.nii.gz'])],
    #       ['MAS021411',dict(loc_nii =[
    #                                   'epi01_mcf.nii.gz',
    #                                   'epi11_mcf.nii.gz'
    #                                   ],
    #                         fix_nii = [
    #                             'epi04_mcf.nii.gz',
    #                             'epi07_mcf.nii.gz',
    #                             'epi10_mcf.nii.gz'
    #                             ],
    #                         left_nii =  [
    #                                      'epi03_mcf.nii.gz',
    #                                      'epi06_mcf.nii.gz',
    #                                      'epi09_mcf.nii.gz'
    #                             ],
    #                         right_nii = [
    #                             'epi02_mcf.nii.gz',
    #                             'epi05_mcf.nii.gz',
    #                             'epi08_mcf.nii.gz'
    #                             ])]],

'''
