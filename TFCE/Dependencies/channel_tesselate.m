function faces = channel_tesselate( vertices )
% CHANNEL_TESSELATE: Tesselate a set of EEG or MEG sensors, for display purpose only.
%
% USAGE:  faces = channel_tesselate( vertices )
%
% INPUT:  
%    - vertices : [Nx3], set of 3D points (MEG or EEG sensors)
% OUTPUT:
%    - faces    : [Mx3], result of the tesselation

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2010 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2009-2010

% === TESSELATE ===
% 2D Projection
[X,Y] = bst_project_2d(vertices(:,1), vertices(:,2), vertices(:,3));
% Compute best fitting sphere
bfs_center = bst_bfs(vertices)';
% Center vertices on BFS center
coordC = bst_bsxfun(@minus, vertices, bfs_center);
% Normalize coordinates
coordC = bst_bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2,2)));
% Tesselation of the sensor array
faces  = convhulln(coordC);
% Get border of the representation
border = convhull(X,Y);
% Keep faces inside the border
iInside = ~(ismember(faces(:,1),border) & ismember(faces(:,2),border)& ismember(faces(:,3),border));
faces   = faces(iInside, :);

% === REMOVE UNNECESSARY TRIANGLES ===
% For instance: the holes for the ears on high-density EEG caps
my_norm = @(v)sqrt(sum(v .^ 2, 2));
% Get coordinates of vertices for each face
vertFacesX = reshape(vertices(reshape(faces,1,[]), 1), size(faces));
vertFacesY = reshape(vertices(reshape(faces,1,[]), 2), size(faces));
vertFacesZ = reshape(vertices(reshape(faces,1,[]), 3), size(faces));
% For each face : compute triangle excentricity and perimeter
triSides = [my_norm([vertFacesX(:,1)-vertFacesX(:,2), vertFacesY(:,1)-vertFacesY(:,2), vertFacesZ(:,1)-vertFacesZ(:,2)]), ...
            my_norm([vertFacesX(:,1)-vertFacesX(:,3), vertFacesY(:,1)-vertFacesY(:,3), vertFacesZ(:,1)-vertFacesZ(:,3)]), ...
            my_norm([vertFacesX(:,2)-vertFacesX(:,3), vertFacesY(:,2)-vertFacesY(:,3), vertFacesZ(:,2)-vertFacesZ(:,3)])];
triPerimeter = sum(triSides, 2);
% Threshold values
thresholdPerim = mean(triPerimeter) + 4.5 * std(triPerimeter);
% Apply threshold
faces(triPerimeter > thresholdPerim, :) = [];



    
    