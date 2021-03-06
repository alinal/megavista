function [img2std] = mrAnatComputeCannonicalXformFromDicomXform(xform, imDim)
% Create a transform from scanner space to a standard space
%
%   [img2std] = ...
%    mrAnatComputeCannonicalXformFromDicomXform(scannerToImXform, imDim)
%
% Given a DICOM-standard scanner-to-image space transform matrix (eg. the
% NIFTI qto_xyz xform), this function will return a 4x4 xform that will
% reorient the volume to a standard axial orientation. (use
% applyCannonicalXform to actually do the reorientation).
%
% The img2std xform will reorient axes so that right-left is along the
% x-axis, anterior-posterior is along the y-axis, superior-inferior is
% along the z-axis, and the leftmost???, anterior-most, superior-most point
% is at 0,0,0 (which, for Analyze/NIFTI, is the lower left-hand corner of
% the last slice).  *** Please check whether RAS has the left most at 0 or
% the rightmost at 0. BW ***
%
% Note that the img2std matrix assumes a PRE-* format and that the n 
% coordinates to be transformed are in a 4xn array. Eg:
%   imgCoords = [0 0 0 1; 0 1 0 1; 1 0 0 1]' 
%   stdCoords = img2std*imgCoords
%
% The reorientation involves only cannonical rotations and mirror flips. 
% Also, note that the reorientation depends on the input xform being a
% proper DICOM xform that converts coords from image space to the DICOM
% standard physical space. Obviously, if this info is wrong, the
% reorientation will be wrong.
%
% Web resources:
%    mrvBrowseSVN('mrAnatComputeCannonicalXformFromDicomXform');
%
% SEE ALSO: 
%   applyCannonicalXform
%   computeCannonicalXformFromIfile
%
% HISTORY:
%   2007.04.19 RFD (bob@white.stanford.edu): wrote it based on
%   computeCannonicalXformFromIfile.
%   


% Compute the scanner-space locations of the image volume corners
imgCorners = [1 1 1 1; imDim(1) 1 1 1; 1 imDim(2) 1 1; imDim(1) imDim(2) 1 1; ...
    1 1 imDim(3) 1; imDim(1) 1 imDim(3) 1; 1 imDim(2) imDim(3) 1; imDim(1) imDim(2) imDim(3) 1];
volRas = xform*imgCorners';
volRas = volRas(1:3,:)';
% The same volRas image corners represented with 0,1:
volXyz = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1];


% Now we need to find the correct rotation & slice reordering to bring 
% volXyz into our standard space. We do this by finding the most right, 
% most anterior, and most superior point (ras), the most left, most 
% anterior, and most superior point (las), etc. for the current volume 
% orientation. Note that the NIFTI convention is that negative values 
% are left, posterior and inferior. The code below finds the correct 
% rotation by measuring the distance from each of the 8 corners to a 
% point in space that is very far to the left, superior and anterior 
% (-1000,1000,1000). Then, we find which of the 8 corners is closest to
% that point. 
d = sqrt((-1000-volRas(:,1)).^2 + (1000-volRas(:,2)).^2 + (1000-volRas(:,3)).^2);
las = find(min(d)==d); las = las(1);
d = sqrt((1000-volRas(:,1)).^2 + (1000-volRas(:,2)).^2 + (1000-volRas(:,3)).^2);
ras = find(min(d)==d); ras = ras(1);
d = sqrt((-1000-volRas(:,1)).^2 + (-1000-volRas(:,2)).^2 + (1000-volRas(:,3)).^2);
lps = find(min(d)==d); lps = lps(1);
d = sqrt((-1000-volRas(:,1)).^2 + (1000-volRas(:,2)).^2 + (-1000-volRas(:,3)).^2);
lai = find(min(d)==d); lai = lai(1);

% Now we have the indices into volRas/volXyz of the 4 anatomical 
% reference points- las, ras, lps and lai. Put them into a 4x4 matrix 
% of homogeneous coordinates.
volCoords = [volXyz(las,:),1; volXyz(lps,:),1; volXyz(lai,:),1; volXyz(ras,:),1;];

% Now we define how we *want* things to be be. That is, the x,y,z location 
% that we'd like for the las, the lps, the lai and the ras (in homogeneous 
% coords). For example:
%    stdCoords = [0,0,0,1; 0,-1,0,1; 0,0,-1,1; 1,0,0,1];
% will map A-P to y axis, L-R to x-axis, and S-I to z-axis with bottom left 
% corner of slice 1 as the most left, most anterior, most inferior point.
% If you want a diferent orientation, you should only need to change this line.
stdCoords = [0,0,0,1; 0,-1,0,1; 0,0,-1,1; 1,0,0,1];

% The following will produce an affine transform matrix that tells us how 
% to transform to our standard space. To use this xform matrix, do: 
% stdCoords = img2std*imgCoords (assuming imgCoords is an 4xn array of n 
% homogeneous coordinates).
img2std = (volCoords \ stdCoords)';

% Fix the translations so that mirror-flips are achieved by -1 rotations.
% This obtuse code relies on the fact that our xform is just 0s 1s and -1s.
% For the rotation part ([1:3],[1:3]), each row should have only one
% nonzero value. If that value is -1, then that denotes a mirror flip. So,
% we set the translations for those dimensions to be imDim rather than 0.
% (Note that we sum across the columns to find the correct imDim value.)
% This way, we get a valid index rather than a negative coord.
img2std(sum(img2std([1:3],[1:3])')<0,4) = imDim(sum(img2std([1:3],[1:3]))<0)';
img2std(sum(img2std([1:3],[1:3])')>0,4) = 0;

% Note that we have constructed this transform matrix so that it will 
% only involve 90, 180 or 270 deg rotations by specifying corresponding 
% points from cannonical locations (the corners of the volume- see stdCoords 
% and volCoords).

return;
