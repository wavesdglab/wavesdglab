% mesh.numVerBnd   % Number of nodes on the boundary
% mesh.numVerInt   % Number of nodes inside the domain
% mesh.listVerBnd  % List of nodes on the boundary     [matrix mesh.numVerBnd x 1]
% mesh.listVerInt  % List of nodes inside the domain   [matrix mesh.numVerInt x 1]

% mesh.numEdg      % Number of edges
% mesh.numEdgBnd   % Number of edges on the boundary
% mesh.numEdgInt   % Number of edges inside the domain
% mesh.mapEdgToVer % List of edges                     [matrix mesh.numEdg    x 2]
% mesh.listEdgBnd  % List of edges on the boundary     [matrix mesh.numEdgBnd    ]
% mesh.listEdgInt  % List of edges inside the domain   [matrix mesh.numEdgInt    ]
% mesh.tagEdg      % Physical tag for edges            [matrix mesh.numEdg       ]
% mesh.tagEdgBnd   % Physical tag for edges on the boundary

% mesh.mapTriToEdg % Connectivity Triangle-to-Edge     [matrix mesh.numTri    x 3]

% mesh.mapTriToTri % Connectivity Triangle-to-Triangle
% mesh.mapTriToFac % Connectivity Triangle-to-Face (LocalEdge)
% mesh.mapEdgToTri % Connectivity Edge-to-Triangle
% mesh.mapEdgToFac % Connectivity Edge-to-Face (LocalEdge)

function mesh = buildMeshConnectivity(mesh)

% -------------------------------------------------------------------------
% List of boundary/interior edges + Connectivity
% -------------------------------------------------------------------------

% Connectivity matrix "boundary-node to boundary-node" (with boundary-TAG ... 0 is forbidden)

connexBndToBnd = sparse(mesh.numVer,mesh.numVer);
for edgBnd = 1:mesh.numEdgBnd
    ver = mesh.mapEdgBndToVer(edgBnd,:);
    connexBndToBnd(ver(1),ver(2)) = mesh.tagEdgBndFile(edgBnd);
    connexBndToBnd(ver(2),ver(1)) = mesh.tagEdgBndFile(edgBnd);
end

% Connectivity all-node to all-node (with edge-ID)
% Lise of edge (general, boundary, interior)

mesh.numEdg      = (3*mesh.numTri + mesh.numEdgBnd)/2;
mesh.numEdgInt   = mesh.numEdg - mesh.numEdgBnd;

mesh.mapTriToEdg = zeros(mesh.numTri,3);
mesh.mapEdgToVer = zeros(mesh.numEdg,2);

mesh.listEdgBnd  = zeros(mesh.numEdgBnd,1);
mesh.listEdgInt  = zeros(mesh.numEdgInt,1);
mesh.tagEdg      = zeros(mesh.numEdg,1);
mesh.tagEdgBnd   = zeros(mesh.numEdgBnd,1);

edg    = 0;
edgBnd = 0;
edgInt = 0;
connexVerToVer = sparse(mesh.numVer,mesh.numVer);
for tri = 1:mesh.numTri
    n1 = mesh.mapTriToVer(tri,:);
    n2 = [n1(2) n1(3) n1(1)]';
    for fac = 1:3
        if (connexVerToVer(n1(fac),n2(fac)) ~= 0)
            mesh.mapTriToEdg(tri,fac) = -connexVerToVer(n1(fac),n2(fac));
        else
            edg = edg + 1;
            connexVerToVer(n1(fac),n2(fac)) = edg;
            connexVerToVer(n2(fac),n1(fac)) = edg;
            mesh.mapTriToEdg(tri,fac) = edg;
            mesh.mapEdgToVer(edg,:) = [n1(fac) n2(fac)];
%             if(n1(fac) < n2(fac))
%                 mesh.mapEdgToVer(edg,:) = [n1(fac) n2(fac)];
%             else
%                 mesh.mapEdgToVer(edg,:) = [n2(fac) n1(fac)];
%             end
            if(connexBndToBnd(n1(fac),n2(fac)) > 0)
                edgBnd = edgBnd+1;
                mesh.listEdgBnd(edgBnd) = edg;
                mesh.tagEdgBnd(edgBnd)  = connexBndToBnd(n1(fac),n2(fac));
                mesh.tagEdg(edg)        = connexBndToBnd(n1(fac),n2(fac));
            else
                edgInt = edgInt+1;
                mesh.listEdgInt(edgInt) = edg;
            end
        end
    end
end
assert(edg    == mesh.numEdg   );
assert(edgBnd == mesh.numEdgBnd);
assert(edgInt == mesh.numEdgInt);

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
    ver = mesh.mapEdgToVer(mesh.listEdgBnd(edg),:);
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