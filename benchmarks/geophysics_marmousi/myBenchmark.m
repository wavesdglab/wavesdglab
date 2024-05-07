function mesh = myBenchmark()

global omega nLambda

% -------------------------------------------------------------------------
% Mesh
% -------------------------------------------------------------------------

myData();
linkMsh = 'benchmarks/geophysics_marmousi/myMesh.msh';
linkGeo = 'benchmarks/geophysics_marmousi/myMesh.geo';
system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -setnumber FREQ ' num2str(omega) ' -setnumber N_LAMBDA ' num2str(nLambda)]);
mesh = readMesh2D(linkMsh);

% -------------------------------------------------------------------------
% Boundary condition
% -------------------------------------------------------------------------

global edgTagToBC
bndTag = {2000};
BC = {'ABC'};
edgTagToBC = containers.Map(bndTag,BC);

% -------------------------------------------------------------------------
% Source point
% -------------------------------------------------------------------------

global pntSouTag
pntSouTag = 1000;

% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

% BP model
% Source:
%   https://www.geoazur.fr/WIND/bin/view/Main/Data/WebHome
% Reference:
%   F. J. Billette and S. Brandsberg-Dahl, The 2004 BP Velocity Benchmark,
%   Extended Abstracts, 67th Annual EAGE Conference & Exhibition, Madrid, Spain, 2004

Ix = 2301;
Iy = 751;
dx = 4; % meter
dy = 4; % meter
Lx = Ix*dx;
Ly = Iy*dy;

fileVelocity = fopen('benchmarks/geophysics_marmousi/data/vp.bin');
fileDensity = fopen('benchmarks/geophysics_marmousi/data/rho.bin');
dataVelocity = fread(fileVelocity,[Iy Ix],'single');
dataDensity = fread(fileDensity,[Iy Ix],'single');

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
    i = floor(x/dx)+1;
    j = floor(-y/dy)+1;
    cArray(tri) = dataVelocity(j,i);
    rhoArray(tri) = dataDensity(j,i);
    
    cArray(tri) = 5000;
    rhoArray(tri) = 1;
    % cArray(tri) = 
    % i0 = ceil(x/dx)+1;
    % i1 = floor(x/dx)+1;
    % j0 = ceil(-y/dy)+1;
    % j1 = floor(-y/dy)+1;
    % cArray(tri) = mean(dataVelocity([j0 j1], [i0 i1]),'all');
    % rhoArray(tri) = mean(dataDensity([j0 j1], [i0 i1]),'all');
    etaArray(tri) = rhoArray(tri) * cArray(tri);
    kArray(tri) = omega / cArray(tri);
end

end