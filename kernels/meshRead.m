function mesh = meshRead(namefile)

% global mesh.numVer      % Number of nodes                                              -- OK
% global mesh.numVerBnd   % Number of nodes on the boundary
% global mesh.numVerInt   % Number of nodes inside the domain
% global mesh.listVerBnd  % List of nodes on the boundary   [matrix mesh.numVerBnd x 1]
% global mesh.listVerInt  % List of nodes inside the domain [matrix mesh.numVerInt x 1]
% global mesh.coord       % Coordinates of vertices         [matrix mesh.numVer    x 2]  -- OK

% global mesh.numEdg      % Number of edges
% global mesh.numEdgBnd   % Number of edges on the boundary
% global mesh.numEdgInt   % Number of edges inside the domain
% global mesh.listEdg     % List of edges                   [matrix mesh.numEdg    x 2]
% global mesh.listEdgBnd  % List of edges on the boundary   [matrix mesh.numEdgBnd    ]
% global mesh.listEdgInt  % List of edges inside the domain [matrix mesh.numEdgInt    ]
% global mesh.tagEdg      % Physical tag for edges          [matrix mesh.numEdg       ]
% global mesh.tagEdgBnd   % Physical tag for edges on the boundary

% global mesh.numTri      % Number of triangles                                          -- OK
% global mesh.mapTriToVer % Connectivity Triangle-to-Vertex     [matrix mesh.numTri x 3] -- OK
% global mesh.mapTriToEdg % Connectivity Triangle-to-Edge       [matrix mesh.numTri x 3]

% global mesh.mapTriToTri % Connectivity Triangle-to-Triangle
% global mesh.mapTriToFac % Connectivity Triangle-to-Face (LocalEdge)
% global mesh.mapEdgToTri % Connectivity Edge-to-Triangle
% global mesh.mapEdgToFac % Connectivity Edge-to-Face (LocalEdge)

% -------------------------------------------------------------------------
% Read the file
% -------------------------------------------------------------------------

file = fopen(namefile,'r');
if (file <= 0)
    error(['Mesh ' namefile ' not found!']);
end

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
NumElements = str2num(fgetl(file));
mesh.mapTriToVer = [];
TagEdgBndFile    = [];
ListEdgBndFile   = [];
for i = 1:NumElements
    line = str2num(fgetl(file));
    % LIN
    if (line(2) == 1)
        TagEdgBndFile  = [TagEdgBndFile; line(4)];
        ListEdgBndFile = [ListEdgBndFile; line(6:7)];
    end
    % TRI
    if (line(2) == 2)
        mesh.mapTriToVer = [mesh.mapTriToVer; line(end-2:end)];
    end
end
NumEdgBndFile = size(ListEdgBndFile,1);
mesh.numTri = size(mesh.mapTriToVer,1);

% Close the file

fclose(file);

% fprintf('---------------------------------------------------------\n');
% fprintf('Mesh    : mesh.numTri    = %i\n', mesh.numTri);

% -------------------------------------------------------------------------
% List of boundary/interior edges + Connectivity
% -------------------------------------------------------------------------

% Connectivity boundary-node to boundary-node (with boundary-TAG ... 0 is forbidden)

mapBndToBnd = sparse(mesh.numVer,mesh.numVer);
for edg = 1:NumEdgBndFile
    ver = ListEdgBndFile(edg,:);
    mapBndToBnd(ver(1),ver(2)) = TagEdgBndFile(edg);
    mapBndToBnd(ver(2),ver(1)) = TagEdgBndFile(edg);
end

% Connectivity all-node to all-node (with edge-ID)
% Lise of edge (general, boundary, interior)

mesh.mapTriToEdg = zeros(mesh.numTri,3);
mesh.numEdg     = 0;
mesh.numEdgBnd  = 0;
mesh.numEdgInt  = 0;
mesh.listEdg    = [];
mesh.listEdgBnd = [];
mesh.listEdgInt = [];
mesh.tagEdgBnd  = [];

edg = 0;
mapVerToVer = sparse(mesh.numVer,mesh.numVer);
for tri = 1:mesh.numTri
    n1 = mesh.mapTriToVer(tri,:);
    n2 = [n1(2) n1(3) n1(1)]';
    for fac = 1:3
        if (mapVerToVer(n1(fac),n2(fac)) ~= 0)
            mesh.mapTriToEdg(tri,fac) = -mapVerToVer(n1(fac),n2(fac));
        else
            mesh.listEdg = [mesh.listEdg; [n1(fac) n2(fac)]];
            edg = edg + 1;
            mapVerToVer(n1(fac),n2(fac)) = edg;
            mapVerToVer(n2(fac),n1(fac)) = edg;
            mesh.mapTriToEdg(tri,fac) = edg;
            if(mapBndToBnd(n1(fac),n2(fac)) > 0)
                mesh.numEdgBnd  = mesh.numEdgBnd+1;
                mesh.listEdgBnd = [mesh.listEdgBnd; edg];
                mesh.tagEdgBnd  = [mesh.tagEdgBnd; mapBndToBnd(n1(fac),n2(fac))];
            else
                mesh.numEdgInt  = mesh.numEdgInt+1;
                mesh.listEdgInt = [mesh.listEdgInt; edg];
            end
        end
    end
end
mesh.numEdg = edg;

mesh.tagEdg = zeros(mesh.numEdg,1);
for edgBnd = 1:NumEdgBndFile
    mesh.tagEdg(mesh.listEdgBnd(edgBnd)) = mesh.tagEdgBnd(edgBnd);
end

% Connectivity edge to triangle/face
% Connectivity triangle to triangle/face

mesh.mapTriToTri = zeros(mesh.numTri,3);
mesh.mapTriToFac = zeros(mesh.numTri,3);
mesh.mapEdgToTri = zeros(mesh.numEdg,2);
mesh.mapEdgToFac = zeros(mesh.numEdg,2);

for tri = 1:mesh.numTri
    for fac = 1:3
        edg = abs(mesh.mapTriToEdg(tri,fac));
        if(mesh.mapEdgToTri(edg,1) == 0)
            mesh.mapEdgToTri(edg,1) = tri;
            mesh.mapEdgToFac(edg,1) = fac;
        else
            mesh.mapEdgToTri(edg,2) = tri;
            mesh.mapEdgToFac(edg,2) = fac;
            tri2 = mesh.mapEdgToTri(edg,1);
            fac2 = mesh.mapEdgToFac(edg,1);
            mesh.mapTriToTri(tri,fac) = tri2;
            mesh.mapTriToFac(tri,fac) = fac2;
            mesh.mapTriToTri(tri2,fac2) = tri;
            mesh.mapTriToFac(tri2,fac2) = fac;
        end
    end
end

% fprintf('Mesh    : mesh.numEdg    = %i\n', mesh.numEdg);
% fprintf('Mesh    : mesh.numEdgBnd = %i\n', mesh.numEdgBnd);
% fprintf('Mesh    : mesh.numEdgInt = %i\n', mesh.numEdgInt);
% fprintf('Mesh    : (3*mesh.numTri = mesh.numEdgBnd + 2*mesh.numEdgInt) -- (%i = %i)\n', 3*mesh.numTri, mesh.numEdgBnd + 2*mesh.numEdgInt);

% -------------------------------------------------------------------------
% List of boundary/interior nodes
% -------------------------------------------------------------------------

maskVerBnd = zeros(mesh.numVer,1);
for edg=1:mesh.numEdgBnd
    ver = mesh.listEdg(mesh.listEdgBnd(edg),:);
    maskVerBnd(ver(1)) = 1;
    maskVerBnd(ver(2)) = 1;
end
mesh.listVerBnd = find(maskVerBnd);
mesh.listVerInt = find(~maskVerBnd);
mesh.numVerBnd = size(mesh.listVerBnd,1);
mesh.numVerInt = size(mesh.listVerInt,1);

% fprintf('Mesh    : mesh.numVer    = %i\n', mesh.numVer);
% fprintf('Mesh    : mesh.numVerBnd = %i\n', mesh.numVerBnd);
% fprintf('Mesh    : mesh.numVerInt = %i\n', mesh.numVerInt);
% fprintf('---------------------------------------------------------\n');

end