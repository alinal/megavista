function mfmSaveStuff(unfoldMesh, params, meshCurvature, gLocs2d, gLocs3d, unfolded2D, numNodes,...
    layerCurvature, startTime, endTime, maps)
% After unfolding a mesh into a flat map, save a bunch of stuff
%
%   mfmSaveStuff(unfoldMesh, params, meshCurvature, gLocs2d, gLocs3d, numNodes,...
%       layerCurvature, startTime, endTime, maps)
%
% Author: Winawer
%  
%   Sub-routine derived from Alex's unfoldMeshFromGUI code.
%
% See Also:  unfoldMeshFromGUI


% Get all the variables we may need
meshFileName     = params.meshFileName;
grayFileName     = params.grayFileName;
flatFileName     = params.flatFileName;
startCoords      = params.startCoords;
scaleFactor      = params.scaleFactor;
perimDist        = params.perimDist;
statusHandle     = params.statusHandle;
busyHandle       = params.busyHandle;
spacingMethod    = params.spacingMethod;
adjustSpacing    = params.adjustSpacing;
gridSpacing      = params.gridSpacing;
showFigures      = params.showFigures;
saveExtra        = params.saveExtra;
truePerimDist    = params.truePerimDist;
hemi             = params.hemi;
nperims          = params.NPERIMS;
saveIntermeidate = params.SAVE_INTERMEDIATE;
numberOfSteps    = params.NUMBEROFSTEPS;

statusStringAdd(statusHandle,'Saving:');
statusStringAdd(statusHandle,flatFileName);

messageString=sprintf('Unfold started at %s\nFinished at %s',datestr(startTime),datestr(endTime));
statusStringAdd(statusHandle,messageString);

% Always show the curvature map, even when the others are suppressed.
mRange=linspace(-perimDist,perimDist,256);
unfoldPlotCurvatureMap(gLocs2d,numNodes,mRange,layerCurvature); 

if (showFigures)
	% save the curvature map
	[p f ext] = fileparts(flatFileName);
	savePath = fullfile(p, [f ' Curvature Map.png']);
	saveas(gcf, savePath);
	fprintf('[%s]: Saved %s.\n', mfilename, savePath);
end

ensureDirExists( fileparts(flatFileName) );

if (saveExtra)
    %try % we added "try and catch" to sweep some error under the rug...
        statusStringAdd(statusHandle,'Saving user data.');
        statusString=char(get(statusHandle,'UserData'));

        % This is the number of l1,l2,l3 and l4 nodes in gLocs3d.
        % Useful for identifying the layers of the
        % various gLocs3d points (since gLocs3d is generated by concatenating l1,l2,l3,l4 gLocs3d)
        infoStr.numNodes.num=numNodes;
        infoStr.numNodes.comment='This is the number of L1,L2,L3 and L4 nodes in gLocs3d. useful for identifying the layers of the various gLocs3d points (since gLocs3d is generated by concatenating L1,L2,L3,L4 gLocs3d)';

        % Save area error maps
        infoStr.faceArea.areaList3D=maps.areaList3D;
        infoStr.faceArea.areaList2D=maps.areaList2D;
        infoStr.faceArea.uniqueFaces=unfoldMesh.uniqueFaceIndexList;
        infoStr.faceArea.errorMap=maps.areaErrorMap;
        infoStr.faceArea.comment='Error map calculated using the center of gravity of the faces and the 3D area/2D area of each face';
        infoStr.faceArea.faceCOGs=[maps.meanY,maps.meanX];
        infoStr.faceArea.errorList=maps.errorList;
        infoStr.faceArea.originalUnfoldMeshVertexList=unfoldMesh.uniqueVertices;

        infoStr.perimDist=perimDist;

        infoStr.startTime=datestr(startTime);
        infoStr.endTime=datestr(endTime);
        infoStr.perimType=truePerimDist;

        infoStr.meshFile=meshFileName;
        infoStr.grayFile=grayFileName;

        unfoldMeshSummary.startCoords=startCoords;
        unfoldMeshSummary.connectionMatrix = unfoldMesh.connectionMatrix;
        unfoldMeshSummary.uniqueVertices = unfoldMesh.uniqueVertices;
        unfoldMeshSummary.uniqueFaceIndexList = unfoldMesh.uniqueFaceIndexList;
        unfoldMeshSummary.internalNodes = unfoldMesh.internalNodes;
        unfoldMeshSummary.orderedUniquePerimeterPoints = unfoldMesh.orderedUniquePerimeterPoints;
        unfoldMeshSummary.scaleFactor = scaleFactor;
        unfoldMeshSummary.locs2d = unfolded2D(:,1:2);
        %unfoldMeshSummary.fullBoundedL1toMeshIndices=fullBoundedL1toMeshIndices; % These are indices into the layer 1 gNodes that each mesh point maps to.

        % we convert the mrGray color index (0-255) to a -1 to 1 curvature value:
        unfoldMeshSummary.curvature = meshCurvature;

        save (flatFileName,'gLocs2d','gLocs3d','meshCurvature','statusString','infoStr','unfoldMeshSummary');
%     catch
%         save (flatFileName,'gLocs2d','gLocs3d','meshCurvature');
%     end

else
    save (flatFileName,'gLocs2d','gLocs3d','meshCurvature');
end

str =sprintf('\n****** End mrFlatMesh  %s *****\n',datestr(endTime));
statusStringAdd(statusHandle,str);

fprintf(1, str);

%----------------------------------
% ---- SUBROUTINES
%----------------------------------
%----------------------------------%----------------------------------
function h = unfoldPlotCurvatureMap(gLocs2d,numNodes,mRange,layerCurvature)

[y x]=meshgrid(mRange);

warning off MATLAB:griddata:DuplicateDataPoints;
fl = griddata(gLocs2d(1:numNodes(1),1),gLocs2d(1:numNodes(1),2),layerCurvature{1},x,y);
warning on MATLAB:griddata:DuplicateDataPoints;

h = figure;  
imagesc((rot90(fl))); colormap gray; title('Curvature map'); axis image

return;