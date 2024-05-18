function mesh = myBenchmark()

% -------------------------------------------------------------------------
% PML
% -------------------------------------------------------------------------

% Boundary condition for the exterior boundary of the PML
global edgTagToBC
edgTag = {201};
BC = {'NEU'};
edgTagToBC = containers.Map(edgTag,BC);

% Perfectly matched layer
global PML_TYPE
if(isempty(PML_TYPE))
    PML_TYPE = 'Circular';
end

% -------------------------------------------------------------------------
% Geometry and Mesh
% -------------------------------------------------------------------------

global h LdomX LdomY LpmlX LpmlY Rdisk Rdom Rpml

switch PML_TYPE
    case 'Rectangular'
        
        if(isempty(LdomX)) LdomX = 1.3; end
        if(isempty(LdomY)) LdomY = 1.3; end
        if(isempty(LpmlX)) LpmlX = 0.2; end
        if(isempty(LpmlY)) LpmlY = 0.2; end
        if(isempty(Rdisk)) Rdisk = 1; end
        
        linkMsh = 'output/mesh.msh';
        linkGeo = 'benchmarks/scattering_disk_penetrable/scattDiskPenetrable_RectangularPML.geo';
        system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ...
            ' -setnumber LdomX ' num2str(LdomX) ' -setnumber LdomY ' num2str(LdomY) ...
            ' -setnumber LpmlX ' num2str(LpmlX) ' -setnumber LpmlY ' num2str(LpmlY) ...
            ' -setnumber Rdisk ' num2str(Rdisk)]);
        mesh = readMesh2D(linkMsh);
        
    case 'Circular'
        
        if(isempty(Rdisk)) Rdisk = 1; end
        if(isempty(Rdom)) Rdom = 1.5; end
        if(isempty(Rpml)) Rpml = 0.5; end
        
        linkMsh = 'output/mesh.msh';
        linkGeo = 'benchmarks/scattering_disk_penetrable/scattDiskPenetrable_CircularPML.geo';
        system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ...
            ' -setnumber Rdisk ' num2str(Rdisk) ...
            ' -setnumber Rdom ' num2str(Rdom) ...
            ' -setnumber Rpml ' num2str(Rpml)]);
        mesh = readMesh2D(linkMsh);
        
end

% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

% Physical parameters for the benchmark
global omega cAir cObj rhoAir rhoObj

% Define tables of coefficients
global rhoArray cArray etaArray kArray
rhoArray = zeros(mesh.numTri,1);
cArray = zeros(mesh.numTri,1);
etaArray = zeros(mesh.numTri,1);
kArray = zeros(mesh.numTri,1);

for tri = 1:mesh.numTri
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    x = (V1(1,1)+V2(1,1)+V3(1,1))/3;
    y = (V1(1,2)+V2(1,2)+V3(1,2))/3;
    if(sqrt(x*x+y*y) < Rdisk)
        cArray(tri) = cObj;
        rhoArray(tri) = rhoObj;
    else
        cArray(tri) = cAir;
        rhoArray(tri) = rhoAir;
    end
    etaArray(tri) = rhoArray(tri) * cArray(tri);
    kArray(tri) = omega / cArray(tri);
end

end