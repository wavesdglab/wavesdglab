function readMesh(namefile)

global NumNod      % Number of nodes
global NumNodBnd   % Number of nodes on the boundary
global NumNodInt   % Number of nodes inside the domain
global ListNodBnd  % List of nodes on the boundary   [matrix NumNodBnd x 1]
global ListNodInt  % List of nodes inside the domain [matrix NumNodInt x 1]
global Coord       % Coordinates des sommets         [matrix NumNod    x 2]

global NumEdg      % Number of edges
global NumEdgBnd   % Number of edges on the boundary
global NumEdgInt   % Number of edges inside the domain
global ListEdg     % List of edges                   [matrix NumEdg    x 2]
global ListEdgBnd  % List of edges on the boundary   [matrix NumEdgBnd    ]
global ListEdgInt  % List of edges inside the domain [matrix NumEdgInt    ]
global TagEdgBnd   % Physical tag for edges on the boundary

global NumTri      % Number of triangles
global ListTri     % List of triangles               [matrix NumTri    x 3]

global mapTriToEdg % Connectivity Triangle-to-Edge
global mapTriToTri % Connectivity Triangle-to-Triangle
global mapTriToFac % Connectivity Triangle-to-Face (LocalEdge)
global mapEdgToTri % Connectivity Edge-to-Triangle
global mapEdgToFac % Connectivity Edge-to-Face (LocalEdge)

% -------------------------------------------------------------------------
% Read the file
% -------------------------------------------------------------------------

file = fopen(namefile,'r');
if (file <= 0)
    error(['Mesh ' namefile ' not found!']);
end

% Read nodes

while (~strcmp(fgetl(file),'$Nodes')) end
NumNod  = str2num(fgetl(file));
Coord   = zeros(NumNod,2);
for i=1:NumNod
    line = str2num(fgetl(file));
    Coord(i,:) = line(2:3);
end

% Read elements

while (~strcmp(fgetl(file),'$Elements')) end
NumElements = str2num(fgetl(file));
ListTri        = [];
TagEdgBndFile  = [];
ListEdgBndFile = [];
for i = 1:NumElements
    line = str2num(fgetl(file));
    if (line(2) == 1)
        TagEdgBndFile  = [TagEdgBndFile; line(4)];
        ListEdgBndFile = [ListEdgBndFile; line(6:7)];
    end
    if (line(2) == 2)
        ListTri = [ListTri; line(end-2:end)];
    end
end
NumEdgBndFile = size(ListEdgBndFile,1);
NumTri = size(ListTri,1);

% Close the file

fclose(file);

fprintf('Mesh    : NumTri    = %i\n', NumTri);

% -------------------------------------------------------------------------
% List of boundary/interior edges + Connectivity
% -------------------------------------------------------------------------

% Connectivity boundary-node to boundary-node (with boundary-TAG)

mapBndToBnd = sparse(NumNod,NumNod);
for edg = 1:NumEdgBndFile
    n = ListEdgBndFile(edg,:);
    mapBndToBnd(n(1),n(2)) = TagEdgBndFile(edg);
    mapBndToBnd(n(2),n(1)) = TagEdgBndFile(edg);
end

% Connectivity all-node to all-node (with edge-ID)
% Lise of edge (general, boundary, interior)

mapTriToEdg = zeros(NumTri,3);
NumEdg = 0;
NumEdgBnd = 0;
NumEdgInt = 0;
ListEdg    = [];
ListEdgBnd = [];
ListEdgInt = [];
TagEdgBnd  = [];

edg = 0;
mapNodToNod = sparse(NumNod,NumNod);
for tri = 1:NumTri
    n1 = ListTri(tri,:);
    n2 = [n1(2) n1(3) n1(1)]';
    for fac = 1:3
        if (mapNodToNod(n1(fac),n2(fac)) ~= 0)
            mapTriToEdg(tri,fac) = -mapNodToNod(n1(fac),n2(fac));
        else
            ListEdg = [ListEdg; [n1(fac) n2(fac)]];
            edg = edg + 1;
            mapNodToNod(n1(fac),n2(fac)) = edg;
            mapNodToNod(n2(fac),n1(fac)) = edg;
            mapTriToEdg(tri,fac) = edg;
            if(mapBndToBnd(n1(fac),n2(fac)) > 0)
                NumEdgBnd  = NumEdgBnd+1;
                ListEdgBnd = [ListEdgBnd; edg];
                TagEdgBnd  = [TagEdgBnd; mapBndToBnd(n1(fac),n2(fac))];
            else
                NumEdgInt  = NumEdgInt+1;
                ListEdgInt = [ListEdgInt; edg];
            end
        end
    end
end
NumEdg = edg;

% Connectivity edge to triangle/face
% Connectivity triangle to triangle/face

mapTriToTri = zeros(NumTri,3);
mapTriToFac = zeros(NumTri,3);
mapEdgToTri = zeros(NumEdg,2);
mapEdgToFac = zeros(NumEdg,2);

for tri = 1:NumTri
    for fac = 1:3
        edg = abs(mapTriToEdg(tri,fac));
        if(mapEdgToTri(edg,1) == 0)
            mapEdgToTri(edg,1) = tri;
            mapEdgToFac(edg,1) = fac;
        else
            mapEdgToTri(edg,2) = tri;
            mapEdgToFac(edg,2) = fac;
            tri2 = mapEdgToTri(edg,1);
            fac2 = mapEdgToFac(edg,1);
            mapTriToTri(tri,fac) = tri2;
            mapTriToFac(tri,fac) = fac2;
            mapTriToTri(tri2,fac2) = tri;
            mapTriToFac(tri2,fac2) = fac;
        end
    end
end

fprintf('Mesh    : NumEdg    = %i\n', NumEdg);
fprintf('Mesh    : NumEdgBnd = %i\n', NumEdgBnd);
fprintf('Mesh    : NumEdgInt = %i\n', NumEdgInt);
fprintf('Mesh    : (3*NumTri = NumEdgBnd + 2*NumEdgInt) -- (%i = %i)\n', 3*NumTri, NumEdgBnd + 2*NumEdgInt);

% -------------------------------------------------------------------------
% List of boundary/interior nodes
% -------------------------------------------------------------------------

maskNodBnd = zeros(NumNod,1);
for edg=1:NumEdgBnd
    n = ListEdg(ListEdgBnd(edg),:);
    maskNodBnd(n(1)) = 1;
    maskNodBnd(n(2)) = 1;
end
ListNodBnd = find(maskNodBnd);
ListNodInt = find(~maskNodBnd);
NumNodBnd = size(ListNodBnd,1);
NumNodInt = size(ListNodInt,1);

fprintf('Mesh    : NumNod    = %i\n', NumNod);
fprintf('Mesh    : NumNodBnd = %i\n', NumNodBnd);
fprintf('Mesh    : NumNodInt = %i\n', NumNodInt);

end