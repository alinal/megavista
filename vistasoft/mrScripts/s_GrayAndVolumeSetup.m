%% s_GrayAndVolumeSetup
%
% Illustrates how to initialize the gray and volume views
%
% You might set up a path file that includes vistadata on your path, and
% you might call it vistaDataPath 
%
% Stanford VISTA


%% Initialize the key variables and data path:
anatDir     = fullfile(mrvDataRootPath,'anatomy','anatomyNIFTI'); % Directory containing anatomies
volAnat     = fullfile(anatDir, 't1.nii.gz'); % Path to volume anatomy
volSegm     = fullfile(anatDir, 't1_class.nii.gz'); % Path to segmentation
nGrayLayers = 5;

% You must tell Matlab where the data directory is.
global HOMEDIR
HOMEDIR = fullfile(mrvDataRootPath,'functional','vwfaloc');

%% Set the volume anatomy path and load the view:
setVAnatomyPath(volAnat);    % Set the volume path
vw_vol = initHiddenVolume(); % Initialize a volume view

%% Grow necessary gray layers from volume and load the view:
buildGrayCoords(vw_vol, [], [], {volSegm}, nGrayLayers); % Use it to grow gray layers
vw_gray = initHiddenGray(); % Initialize a gray view

%% END