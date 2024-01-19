function mesh = myBenchmark(h)

global BCObstacle L L_PML BCPML

BCPML = 'DIR';
BCObstacle = 'NEU';

linkMsh = 'benchmarks/scattering_square/square_cavity.msh';
linkGeo = 'benchmarks/scattering_square/square_cavity.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber L ' num2str(L) ' -setnumber L_PML ' num2str(L_PML)]);
mesh = readMesh2D(linkMsh);

end