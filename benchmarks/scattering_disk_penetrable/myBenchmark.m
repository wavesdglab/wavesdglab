function mesh = myBenchmark()

% -------------------------------------------------------------------------
% PML
% -------------------------------------------------------------------------

% Boundary condition for the exterior boundary of the PML
global edgTagToBC
edgTag = {201};
BC = {'NEU0'};
edgTagToBC = containers.Map(edgTag,BC);

% Perfectly matched layer
global PML_TYPE
if(isempty(PML_TYPE))
    PML_TYPE = 'Circular';
end

% -------------------------------------------------------------------------
% Geometry and Mesh
% -------------------------------------------------------------------------

global h LdomX LdomY LpmlX LpmlY Rdisk Rdom Rpml cObj

switch PML_TYPE
    case 'Rectangular'
        
        if(isempty(LdomX)) LdomX = 1.3; end
        if(isempty(LdomY)) LdomY = 1.3; end
        if(isempty(LpmlX)) LpmlX = 0.2; end
        if(isempty(LpmlY)) LpmlY = 0.2; end
        if(isempty(Rdisk)) Rdisk = 1; end

        % Compute finer mesh size for the disk
        hDisk = h * cObj;

        linkMsh = 'output/mesh.msh';
        linkGeo = append(fileparts(mfilename('fullpath')),'/scattDiskPenetrable_RectangularPML.geo');
        system(['gmsh -2 ' linkGeo ' -o ' linkMsh ...
            ' -setnumber LdomX ' num2str(LdomX) ' -setnumber LdomY ' num2str(LdomY) ...
            ' -setnumber LpmlX ' num2str(LpmlX) ' -setnumber LpmlY ' num2str(LpmlY) ...
            ' -setnumber Rdisk ' num2str(Rdisk) ...
            ' -setnumber h ' num2str(h) ...
            ' -setnumber hDisk ' num2str(hDisk)]);
        mesh = readMesh2D(linkMsh);

    case 'Circular'

        if(isempty(Rdisk)) Rdisk = 1; end
        if(isempty(Rdom)) Rdom = 1.2; end
        if(isempty(Rpml)) Rpml = 0.2; end

        % Compute finer mesh size for the disk
        hDisk = h * cObj;

        linkMsh = 'output/mesh.msh';
        linkGeo = append(fileparts(mfilename('fullpath')),'/scattDiskPenetrable_CircularPML.geo');
        system(['gmsh -2 ' linkGeo ' -o ' linkMsh ...
            ' -setnumber Rdisk ' num2str(Rdisk) ...
            ' -setnumber Rdom ' num2str(Rdom) ...
            ' -setnumber Rpml ' num2str(Rpml) ...
            ' -setnumber h ' num2str(h) ...
            ' -setnumber hDisk ' num2str(hDisk)]);
        mesh = readMesh2D(linkMsh);
        
end

% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

% Physical parameters for the benchmark
global omega cAir rhoAir rhoObj

% Define tables of coefficients
global rho c eta k
rho = zeros(mesh.numTri,1);
c = zeros(mesh.numTri,1);
eta = zeros(mesh.numTri,1);
k = zeros(mesh.numTri,1);

for tri = 1:mesh.numTri
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    x = (V1(1,1)+V2(1,1)+V3(1,1))/3;
    y = (V1(1,2)+V2(1,2)+V3(1,2))/3;
    if(sqrt(x*x+y*y) < Rdisk)
        c(tri) = cObj;
        rho(tri) = rhoObj;
    else
        c(tri) = cAir;
        rho(tri) = rhoAir;
    end
    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);
end

end