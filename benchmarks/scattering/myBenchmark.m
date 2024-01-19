function mesh = myBenchmark(h)

global BCObstacle L L_PML R_disk BCPML

BCPML = 'DIR';
BCObstacle = 'NEU';

linkMsh = 'benchmarks/scattering/scatteringPML.msh';
linkGeo = 'benchmarks/scattering/scatteringPML.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber L ' num2str(L) ' -setnumber L_PML ' num2str(L_PML) ' -setnumber R_disk ' num2str(R_disk)]);
mesh = readMesh2D(linkMsh);

end