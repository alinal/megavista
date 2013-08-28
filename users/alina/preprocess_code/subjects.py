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

#rois=['L_V1_0.25', 'L_V2V_0.25', 'L_V3V_0.25', 'L_V4_0.25',  'l_mFus_p3', 'l_PPA_p4',
#      'L_V2D_0.25',  'L_V3D_0.25', 'L_V3A_0.25',  'L_IPS0_0.25', 'L_IPS1_0.25', 'L_IPS2_0.25' ]

#rois=[ 'R_V1_0.25', 'R_V2V_0.25', 'R_V3V_0.25', 'R_V4_0.25', 'r_pFus_p3', 'r_mFus_p3', 'r_PPA_p4',
#      'R_V2D_0.25',  'R_V3D_0.25', 'R_V3A_0.25',  'R_IPS0_0.25', 'R_IPS1_0.25', 'R_IPS2_0.25',
#      'L_V1_0.25', 'L_V2V_0.25', 'L_V3V_0.25', 'L_V4_0.25',  'l_mFus_p3', 'l_PPA_p4',
#      'L_V2D_0.25',  'L_V3D_0.25', 'L_V3A_0.25',  'L_IPS0_0.25', 'L_IPS1_0.25', 'L_IPS2_0.25']


rois=[ 'R_V1_0.25', 'R_V2V_0.25', 'R_V3V_0.25', 'R_V4_0.25', 'R_LO1_0.25', 'R_LO2_0.25',
       'r_LOf_p3_0.25', 'r_pFus_p3_0.25',  'r_mFus_p3_0.25', 'r_PPA_p4_0.25',
      'R_V2D_0.25',  'R_V3D_0.25', 'R_V3A_0.25', 'R_IPS0_0.25', 'R_IPS1_0.25', 'R_IPS2_0.25',
      'R_IPS3_0.25', 'R_IPS4_0.25', 'R_IPS5_0.25',
      'L_V1_0.25', 'L_V2V_0.25', 'L_V3V_0.25', 'L_V4_0.25', 'L_LO1_0.25','L_LO2_0.25',
      'l_pFus_p3_0.25','l_mFus_p3_0.25',  'l_pSTS_p3_0.25', 'l_PPA_p4_0.25',
      'L_V2D_0.25',  'L_V3D_0.25', 'L_V3A_0.25', 'L_IPS0_0.25', 'L_IPS1_0.25', 'L_IPS2_0.25',
      'L_IPS3_0.25', 'L_IPS4_0.25', 'L_IPS5_0.25']

nuisReg=['mcf1_bp.1D', 'mcf2_bp.1D', 'mcf3_bp.1D', 'mcf4_bp.1D','mcf5_bp.1D', 'mcf6_bp.1D', 'lVent_4mm_stc_avgFt.1D',
        'rVent_4mm_stc_avgFt.1D', 'rWM_6mm_stc_avgFt.1D', 'lWM_6mm_stc_avgFt.1D', 'gray_al_stc_avgFt.1D', 'mcf1_bp_deriv.1D', 'mcf2_bp_deriv.1D',
         'mcf3_bp_deriv.1D', 'mcf4_bp_deriv.1D','mcf5_bp_deriv.1D', 'mcf6_bp_deriv.1D','lVent_4mm_stc_avgFt_deriv.1D',
        'rVent_4mm_stc_avgFt_deriv.1D', 'rWM_6mm_stc_avgFt_deriv.1D', 'lWM_6mm_stc_avgFt_deriv.1D', 'gray_al_stc_avgFt_deriv.1D']

subjects = {  'DCA':[['DCA042511',dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             fix_nii = ['epi04_mcf.nii.gz',
                                        'epi07_mcf.nii.gz',
                                        'epi10_mcf.nii.gz'],
                             left_nii =  ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'],
                             right_nii = ['epi02_mcf.nii.gz',
                                          'epi05_mcf.nii.gz',
                                          'epi08_mcf.nii.gz'])],
                ['DCA041111', dict(loc_nii =['epi01_mcf.nii.gz',
                                             'epi11_mcf.nii.gz'],
                                   fix_nii = ['epi04_mcf.nii.gz',
                                              'epi07_mcf.nii.gz',
                                              'epi10_mcf.nii.gz'],
                                   left_nii = ['epi03_mcf.nii.gz',
                                                'epi06_mcf.nii.gz',
                                                'epi09_mcf.nii.gz'],
                                   right_nii = ['epi02_mcf.nii.gz',
                                                'epi05_mcf.nii.gz',
                                                'epi08_mcf.nii.gz'])]]}


'''
    'CG':[['CG011611', dict(loc_nii =['epi01_mcf_stc.nii.gz',
                                      'epi10_mcf_stc.nii.gz'],
                            fix_nii = ['epi04_mcf_stc.nii.gz',
                                       'epi07_mcf_stc.nii.gz'],
                            right_nii =  ['epi03_mcf_stc.nii.gz',
                                          'epi06_mcf_stc.nii.gz',
                                          'epi09_mcf_stc.nii.gz'],
                            left_nii = ['epi02_mcf_stc.nii.gz',
                                        'epi05_mcf_stc.nii.gz',
                                        'epi08_mcf_stc.nii.gz'])],
          ['CG020611', dict(loc_nii = ['epi01_mcf_stc.nii.gz',
                                       'epi11_mcf_stc.nii.gz'],
                            fix_nii = ['epi04_mcf_stc.nii.gz',
                                       'epi07_mcf_stc.nii.gz',
                                       'epi10_mcf_stc.nii.gz'],
                            right_nii = ['epi03_mcf_stc.nii.gz',
                                         'epi06_mcf_stc.nii.gz',
                                         'epi09_mcf_stc.nii.gz'],
                            left_nii =  ['epi02_mcf_stc.nii.gz',
                                         'epi05_mcf_stc.nii.gz',
                                         'epi08_mcf_stc.nii.gz'])]]}

'CHT':[['CHT042111', dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             right_nii = ['epi04_mcf.nii.gz',
                                        'epi07_mcf.nii.gz',
                                        'epi10_mcf.nii.gz'],
                             left_nii =  ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'],
                             fix_nii = ['epi02_mcf.nii.gz',
                                        'epi05_mcf.nii.gz',
                                        'epi08_mcf.nii.gz'])],

            ['CHT051911',dict(loc_nii =['epi01_mcf.nii.gz',
                                        'epi11_mcf.nii.gz'],
                             right_nii = ['epi04_mcf.nii.gz',
                                          'epi07_mcf.nii.gz',
                                          'epi10_mcf.nii.gz'],
                             left_nii = ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'],
                             fix_nii = ['epi02_mcf.nii.gz',
                                        'epi05_mcf.nii.gz',
                                        'epi08_mcf.nii.gz'])]]}



  '''

         #            'SS':[['SS012311', dict(loc_nii =['epi01_mcf.nii.gz',
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
     #       ['SS011011',dict(loc_nii = ['epi01_mcf.nii.gz',
     #                                   'epi10_mcf.nii.gz'],
     #                        fix_nii = ['epi04_mcf.nii.gz',
     #                                   'epi07_mcf.nii.gz'],
     #                        left_nii = ['epi03_mcf.nii.gz',
     #                                    'epi06_mcf.nii.gz',
     #                                    'epi09_mcf.nii.gz'],
     #                        right_nii =  ['epi02_mcf.nii.gz',
     #                                      'epi05_mcf.nii.gz',
     #                                      'epi08_mcf.nii.gz'])]]
     #                                    }

'''
     'SS':[
         ['SS012311', dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             fix_nii = ['epi04_mcf.nii.gz',
                                        'epi07_mcf.nii.gz',
                                        'epi10_mcf.nii.gz'],
                             left_nii =  ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'],
                             right_nii = ['epi02_mcf.nii.gz',
                                          'epi05_mcf.nii.gz',
                                          'epi08_mcf.nii.gz'])],
           ['SS011011',dict(loc_nii = ['epi01_mcf.nii.gz',
                                       'epi10_mcf.nii.gz'],
                            fix_nii = ['epi04_mcf.nii.gz',
                                       'epi07_mcf.nii.gz'],
                            left_nii = ['epi03_mcf.nii.gz',
                                        'epi06_mcf.nii.gz',
                                        'epi09_mcf.nii.gz'],
                            right_nii =  ['epi02_mcf.nii.gz',
                                          'epi05_mcf.nii.gz',
                                          'epi08_mcf.nii.gz'])]],

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

     'WC':[['WC031911', dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             left_nii = ['epi02_mcf.nii.gz',
                                        'epi05_mcf.nii.gz',
                                        'epi08_mcf.nii.gz'],
                             fix_nii =  ['epi04_mcf.nii.gz',
                                          'epi07_mcf.nii.gz',
                                          'epi10_mcf.nii.gz'],
                             right_nii = ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'])],
           ['WC030311', dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             left_nii = ['epi02_mcf.nii.gz',
                                        'epi05_mcf.nii.gz',
                                        'epi08_mcf.nii.gz'],
                             fix_nii =  ['epi04_mcf.nii.gz',
                                          'epi07_mcf.nii.gz',
                                          'epi10_mcf.nii.gz'],
                             right_nii = ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'])]],


           'DCA':[['DCA042511',dict(loc_nii =['epi01_mcf.nii.gz',
                                       'epi11_mcf.nii.gz'],
                             fix_nii = ['epi04_mcf.nii.gz',
                                        'epi07_mcf.nii.gz',
                                        'epi10_mcf.nii.gz'],
                             left_nii =  ['epi03_mcf.nii.gz',
                                          'epi06_mcf.nii.gz',
                                          'epi09_mcf.nii.gz'],
                             right_nii = ['epi02_mcf.nii.gz',
                                          'epi05_mcf.nii.gz',
                                          'epi08_mcf.nii.gz'])],
                ['DCA041111', dict(loc_nii =['epi01_mcf.nii.gz',
                                             'epi11_mcf.nii.gz'],
                                   fix_nii = ['epi04_mcf.nii.gz',
                                              'epi07_mcf.nii.gz',
                                              'epi10_mcf.nii.gz'],
                                   left_nii = ['epi03_mcf.nii.gz',
                                                'epi06_mcf.nii.gz',
                                                'epi09_mcf.nii.gz'],
                                   right_nii = ['epi02_mcf.nii.gz',
                                                'epi05_mcf.nii.gz',
                                                'epi08_mcf.nii.gz'])]],



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
