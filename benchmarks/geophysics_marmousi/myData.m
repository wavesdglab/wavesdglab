function myData()

% Marmousi model
% Source:
%   https://www.geoazur.fr/WIND/bin/view/Main/Data/WebHome#mjx-eqn%3AMarmousi_(2D_-_acoustic_-_isotropic)

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

figure(1);
imagesc(dataVelocity);
title('Velocity');
colorbar;

figure(2);
imagesc(dataDensity);
title('Density');
colorbar;

% Data for structured GMSH data file:
%   Ox Oy Oz
%   Dx Dy Dz
%   nx ny nz

padding = 10;

data = round(dataVelocity(1:padding:end, 1:padding:end)');
fileGmsh = fopen('benchmarks/geophysics_marmousi/_velocityGmsh.txt','w');
fprintf(fileGmsh,'%f %f %f\n', 0, 0, 0);
fprintf(fileGmsh,'%f %f %f\n', Lx/size(data,1), -Ly/size(data,2), 1);
fprintf(fileGmsh,'%i %i %i\n', size(data,1), size(data,2), 1);
for i=1:size(data,1)
    for j=1:size(data,2)
        fprintf(fileGmsh,'%.0f ', data(i,j));
    end
    fprintf(fileGmsh,'\n');
end
fclose(fileGmsh);

end