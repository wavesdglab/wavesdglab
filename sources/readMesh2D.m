% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% mesh.numVer          % Number of nodes
% mesh.coord           % Coordinates of vertices               [matrix mesh.numVer x 2]

% mesh.numTri          % Number of triangles
% mesh.numEdgBnd       % Number of boundary edges
% mesh.mapTriToVer     % Connectivity Triangle-to-Vertices     [matrix mesh.numTri x 3]
% mesh.mapEdgBndToVer  % Connectivity BoundaryEdge-to-Vertices [matrix mesh.numEdgBnd x 2]
% mesh.tagEdgBndFile

% mesh.hmin            % Min length of edge
% mesh.hmax            % Max length of edge

% -------------------------------------------------------------------------
% Mesh reader - General function
% -------------------------------------------------------------------------

function mesh = readMesh2D(namefile)

% Open file
file = fopen(namefile,'r');
if (file <= 0)
    error(['Mesh ' namefile ' not found!']);
end

% Read mesh format
while (~strcmp(fgetl(file),'$MeshFormat'))
end
meshFormat = str2num(fgetl(file));
if(meshFormat(1) ~= 4.1)
    error(['Mesh format ' meshFormat(1) ' not known!']);
end

% Read entities
while (~strcmp(fgetl(file),'$Entities'))
end
data = str2num(fgetl(file));
geo.numPoints = data(1);
geo.numCurves = data(2);
geo.numSurfaces = data(3);
geo.numVolumes = data(4);
for i=1:geo.numPoints
    data = str2num(fgetl(file));
    geo.pointGeoTag(i) = data(1);
    if(data(5) > 0)
        geo.pointPhyTag(i) = data(6);
    else
        geo.pointPhyTag(i) = -1;
    end
end
for i=1:geo.numCurves
    data = str2num(fgetl(file));
    geo.curveGeoTag(i) = data(1);
    geo.curvePhyTag(i) = data(9);
end
for i=1:geo.numSurfaces
    data = str2num(fgetl(file));
    geo.surfaceGeoTag(i) = data(1);
    geo.surfacePhyTag(i) = data(9);
end
for i=1:geo.numVolumes
    data = str2num(fgetl(file));
    geo.volumeGeoTag(i) = data(1);
    geo.volumePhyTag(i) = data(9);
end

% Read nodes
while (~strcmp(fgetl(file),'$Nodes'))
end
data = str2num(fgetl(file));
numEntities = data(1);
mesh.numVer = data(2);
%mesh.nodeTag = zeros(mesh.numVer,1);
%mesh.nodeDim = zeros(mesh.numVer,1);
mesh.coord = zeros(mesh.numVer,2);
n = 0;
for ent=1:numEntities
    data = str2num(fgetl(file));
    entityDim = data(1);
    numNodesInBlock = data(4);
    for i=1:numNodesInBlock
        %mesh.nodeDim(n+i) = entityDim;
    end
    for i=1:numNodesInBlock
        data = str2num(fgetl(file));
        %mesh.nodeTag(n+i) = str2num(fgetl(file));
    end
    for i=1:numNodesInBlock
        data = str2num(fgetl(file));
        mesh.coord(n+i,:) = data(:,1:2);
    end
    n = n+numNodesInBlock;
end

% Read elements
while (~strcmp(fgetl(file),'$Elements'))
end
data = str2num(fgetl(file));
numEntityBlocks = data(1);
numElementsPNT      = 0;
numElementsLIN      = 0;
numElementsTRI      = 0;
mesh.mapTriToVer    = [];
mesh.tagPntFile  = [];
mesh.mapPntToVer = [];
mesh.tagEdgBndFile  = [];
mesh.mapEdgBndToVer = [];
mesh.tagTriFile = [];
for ent=1:numEntityBlocks
    data = str2num(fgetl(file));
    entityDim = data(1);
    entityTag = data(2);
    entityType = data(3);
    numElementsInBlock = data(4);
    if(entityType == 15)
        for i=1:numElementsInBlock
            data = str2num(fgetl(file));
            phyTag = geo.pointPhyTag(geo.pointGeoTag == entityTag);
            mesh.tagPntFile  = [mesh.tagPntFile; phyTag];
            mesh.mapPntToVer = [mesh.mapPntToVer; data(2)];
        end
        numElementsPNT = numElementsPNT + numElementsInBlock;
    end
    if(entityType == 1)
        for i=1:numElementsInBlock
            data = str2num(fgetl(file));
            phyTag = geo.curvePhyTag(geo.curveGeoTag == entityTag);
            mesh.tagEdgBndFile  = [mesh.tagEdgBndFile; phyTag];
            mesh.mapEdgBndToVer = [mesh.mapEdgBndToVer; data(2:3)];
        end
        numElementsLIN = numElementsLIN + numElementsInBlock;
    end
    if(entityType == 2)
        for i=1:numElementsInBlock
            data = str2num(fgetl(file));
            mesh.tagTriFile = [mesh.tagTriFile; data(1)];
            mesh.mapTriToVer = [mesh.mapTriToVer; data(2:4)];
        end
        numElementsTRI = numElementsTRI + numElementsInBlock;
    end
end

mesh.numTri = size(mesh.mapTriToVer,1);
mesh.numEdgBnd = size(mesh.mapEdgBndToVer,1);

% Close file
fclose(file);

% Mesh sizes
mesh.hmin = 1e10;
mesh.hmax = 0;
for i = 1:mesh.numTri
    V = mesh.mapTriToVer(i,:);
    X = mesh.coord(V,:);
    l1 = norm(X(1,:)-X(2,:));
    l2 = norm(X(2,:)-X(3,:));
    l3 = norm(X(3,:)-X(1,:));
    mesh.hmax = max(max(max(mesh.hmax,l1),l2),l3);
    mesh.hmin = min(min(min(mesh.hmin,l1),l2),l3);
end

end