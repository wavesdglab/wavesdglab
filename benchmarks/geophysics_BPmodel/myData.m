function myData()

% BP model
% Source:
%   https://www.geoazur.fr/WIND/bin/view/Main/Data/WebHome#mjx-eqn%3A2004_BP_salt_model_(2D_-_acoustic_-_isotropic)

Ix = 5395;
Iy = 1911;
dx = 12.5; % meter
dy = 6.25; % meter
Lx = Ix*dx; % 67437.5
Ly = Iy*dy; % 11943.75

fileVelocity = fopen('benchmarks/geophysics_BPmodel/data/vp.bin');
dataVelocity = fread(fileVelocity,[Iy Ix],'single');

% fileDensity = fopen('benchmarks/geophysics_BPmodel/data/rho.bin');
% dataDensity = fread(fileDensity,[Iy Ix],'single');

% figure(1);
% imagesc(dataVelocity);
% title('Velocity');
% colorbar;
% 
% figure(2);
% imagesc(dataDensity);
% title('Density');
% colorbar;

% Data for structured GMSH data file:
%   Ox Oy Oz
%   Dx Dy Dz
%   nx ny nz

padding = 10;

data = round(dataVelocity(1:padding:end, 1:padding:end)');
fileGmsh = fopen('output/velocityGmsh.txt','w');
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