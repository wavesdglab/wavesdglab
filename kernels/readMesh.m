% global mesh.numVer      % Number of nodes                                                -- OK
% global mesh.numVerBnd   % Number of nodes on the boundary
% global mesh.numVerInt   % Number of nodes inside the domain
% global mesh.listVerBnd  % List of nodes on the boundary     [matrix mesh.numVerBnd x 1]
% global mesh.listVerInt  % List of nodes inside the domain   [matrix mesh.numVerInt x 1]
% global mesh.coord       % Coordinates of vertices           [matrix mesh.numVer    x 2]  -- OK

% global mesh.numEdg      % Number of edges
% global mesh.numEdgBnd   % Number of edges on the boundary
% global mesh.numEdgInt   % Number of edges inside the domain
% global mesh.listEdg     % List of edges                     [matrix mesh.numEdg    x 2]
% global mesh.listEdgBnd  % List of edges on the boundary     [matrix mesh.numEdgBnd    ]
% global mesh.listEdgInt  % List of edges inside the domain   [matrix mesh.numEdgInt    ]
% global mesh.tagEdg      % Physical tag for edges            [matrix mesh.numEdg       ]
% global mesh.tagEdgBnd   % Physical tag for edges on the boundary

% global mesh.numTri      % Number of triangles                                           -- OK
% global mesh.mapTriToVer % Connectivity Triangle-to-Vertex   [matrix mesh.numTri    x 3] -- OK
% global mesh.mapTriToEdg % Connectivity Triangle-to-Edge     [matrix mesh.numTri    x 3]

% global mesh.mapTriToTri % Connectivity Triangle-to-Triangle
% global mesh.mapTriToFac % Connectivity Triangle-to-Face (LocalEdge)
% global mesh.mapEdgToTri % Connectivity Edge-to-Triangle
% global mesh.mapEdgToFac % Connectivity Edge-to-Face (LocalEdge)

% =========================================================================
% Mesh reader - General function
% =========================================================================

function mesh = readMesh(namefile)
disp(['--- CALL readMesh']);

file = fopen(namefile,'r');
if (file <= 0)
    error(['Mesh ' namefile ' not found!']);
end

while (~strcmp(fgetl(file),'$MeshFormat'))
end
meshFormat = str2num(fgetl(file));

switch meshFormat(1)
    case 2.2
        mesh = readMeshGmsh2v2(file);
    case 4.1
        mesh = readMeshGmsh4v1(file);
    otherwise
        error(['Mesh format ' meshFormat(1) ' not known!']);
end

fclose(file);

end

% =========================================================================
% Mesh reader - Gmsh format 2.2
% =========================================================================

function mesh = readMeshGmsh2v2(file)
disp(['--- CALL readMeshGmsh2v2']);

% Read nodes

while (~strcmp(fgetl(file),'$Nodes')) end
mesh.numVer = str2num(fgetl(file));
mesh.coord  = zeros(mesh.numVer,2);
for i=1:mesh.numVer
    line = str2num(fgetl(file));
    mesh.coord(i,:) = line(2:3);
end

% Read elements

while (~strcmp(fgetl(file),'$Elements')) end
numElements = str2num(fgetl(file));
mesh.mapTriToVer    = [];
mesh.tagEdgBndFile  = [];
mesh.listEdgBndFile = [];
for i = 1:numElements
    line = str2num(fgetl(file));
    % LIN
    if (line(2) == 1)
        mesh.tagEdgBndFile  = [mesh.tagEdgBndFile;  line(4)];
        mesh.listEdgBndFile = [mesh.listEdgBndFile; line(6:7)];
    end
    % TRI
    if (line(2) == 2)
        mesh.mapTriToVer = [mesh.mapTriToVer; line(end-2:end)];
    end
end
mesh.numTri = size(mesh.mapTriToVer,1);

end

% =========================================================================
% Mesh reader - Gmsh format 4.1
% =========================================================================

function mesh = readMeshGmsh4v1(file)
disp(['--- CALL readMeshGmsh4v1']);

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
    geo.pointTag(i) = data(1);
    geo.X(i) = data(2);
    geo.Y(i) = data(3);
    geo.Z(i) = data(4);
end
for i=1:geo.numCurves
    data = str2num(fgetl(file));
    geo.curveGeoTag(i) = data(1);
    geo.curveGeoTag(i) = data(8);
end
for i=1:geo.numSurfaces
    data = str2num(fgetl(file));
    geo.surfaceGeoTag(i) = data(1);
    geo.surfacePhyTag(i) = data(8);
end
for i=1:geo.numVolumes
    data = str2num(fgetl(file));
    geo.volumeGeoTag(i) = data(1);
    geo.volumePhyTag(i) = data(8);
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
numElementsLIN     = 0;
numElementsTRI      = 0;
mesh.mapTriToVer    = [];
mesh.tagEdgBndFile  = [];
mesh.listEdgBndFile = [];
for ent=1:numEntityBlocks
    data = str2num(fgetl(file));
    entityDim = data(1);
    entityTag = data(2);
    entityType = data(3);
    numElementsInBlock = data(4);
    if(entityType == 1)
        for i=1:numElementsInBlock
            data = str2num(fgetl(file));
            mesh.tagEdgBndFile  = [mesh.tagEdgBndFile; entityTag];
            mesh.listEdgBndFile = [mesh.listEdgBndFile; data(2:3)];
        end
        numElementsLIN = numElementsLIN + numElementsInBlock;
    end
    if(entityType == 2)
        for i=1:numElementsInBlock
            data = str2num(fgetl(file));
            mesh.mapTriToVer = [mesh.mapTriToVer; data(2:4)];
        end
        numElementsTRI = numElementsTRI + numElementsInBlock;
    end
end

mesh.numTri = size(mesh.mapTriToVer,1);

end