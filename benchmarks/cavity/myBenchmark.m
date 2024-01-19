function mesh = myBenchmark(h)

global BCWest BCNorth BCEast BCSouth

BCWest  = 'DIR';
BCNorth = 'DIR';
BCEast  = 'DIR';
BCSouth = 'DIR';

linkMsh = 'benchmarks/cavity/cavity.msh';
linkGeo = 'benchmarks/cavity/cavity.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end