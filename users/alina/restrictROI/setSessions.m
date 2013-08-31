% Get sessions


%sessions={'CG011611' 'CG020611' 'DCA041111' 'DCA042511' ...
%		'SS011011' 'SS012311' 'CHT042111' 'CHT051911'...	
%       'WC030311' 'WC031911'}; 
    
sessions={'SS011011' 'SS012311'};
    
rois={{'r_IOG_p3.mat', 'r_LOf_p3.mat', 'r_pFus_p3.mat', 'r_mFus_p3.mat', 'r_PPA_p4.mat', 'r_pSTS_p3.mat', 'l_IOG_p3.mat',... 
    'l_LOf_p3.mat', 'l_pFus_p3.mat', 'l_mFus_p3.mat', 'l_PPA_p4.mat', 'l_pSTS_p3.mat'}...
    {'r_IOG_p3.mat', 'r_LOf_p3.mat', 'r_pFus_p3.mat', 'r_mFus_p3.mat', 'r_PPA_p4.mat', 'r_pSTS_p3.mat', 'l_IOG_p3.mat',... 
    'l_LOf_p3.mat', 'l_pFus_p3.mat', 'l_mFus_p3.mat', 'l_PPA_p4.mat', 'l_pSTS_p3.mat'}}

% rois={
% {'L_LGN_0.3.mat' 'R_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_V4.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3D.mat' 'L_V3V.mat' 'L_V3A.mat' 'L_VO1.mat' ...
% 'R_IPS0.mat' 'R_IPS1.mat' 'R_IPS2.mat' 'R_LO1.mat' 'R_LO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V4.mat' 'R_V3A.mat' 'R_VO1.mat'}... 
% {'L_LGN_0.3.mat' 'R_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_V4.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3D.mat' 'L_V3V.mat' 'L_V3A.mat' 'L_VO1.mat' ...
% 'R_IPS0.mat' 'R_IPS1.mat' 'R_IPS2.mat' 'R_LO1.mat' 'R_LO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V4.mat' 'R_V3A.mat' 'R_VO1.mat'} ... 
% {'L_LGN_0.3.mat' 'R_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_SPL1.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_SPL1.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}...
% {'L_LGN_0.3.mat' 'R_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_SPL1.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_SPL1.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat'... 
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3D.mat' 'L_V3V.mat' 'L_V3A.mat' ...
% 'R_IPS0.mat' 'R_IPS1.mat' 'R_IPS2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat' } ...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat'... 
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3D.mat' 'L_V3V.mat' 'L_V3A.mat' ...
% 'R_IPS0.mat' 'R_IPS1.mat' 'R_IPS2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat' } ...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_SPL1.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_SPL1.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_SPL1.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_SPL1.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}...
% {'R_LGN_0.3.mat' 'L_LGN_0.3.mat'...
% 'L_IPS0.mat', 'L_IPS1.mat' 'L_IPS2.mat' 'L_IPS3.mat' 'L_IPS4.mat'...
% 'L_IPS5.mat' 'L_LO1.mat' 'L_LO2.mat' 'L_TO1.mat' 'L_TO2.mat'...
% 'L_V1.mat' 'L_V2D.mat' 'L_V2V.mat' 'L_V3A.mat' 'L_V3D.mat' 'L_V3V.mat' 'R_IPS0.mat'...
% 'R_IPS1.mat' 'R_IPS2.mat' 'R_IPS3.mat' 'R_IPS4.mat' 'R_IPS5.mat' 'R_LO1.mat'...
% 'R_LO2.mat' 'R_TO1.mat' 'R_TO2.mat' 'R_V1.mat' 'R_V2D.mat'...
% 'R_V2V.mat' 'R_V3D.mat' 'R_V3V.mat' 'R_V3A.mat'}};
% 
% 
% 
% 
