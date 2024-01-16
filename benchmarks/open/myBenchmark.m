function mesh = myBenchmark(h)

global BCWest BCNorth BCEast BCSouth

BCWest  = 'ABC';
BCNorth = 'ABC';
BCEast  = 'ABC';
BCSouth = 'ABC';

linkMsh = 'benchmarks/open/open.msh';
linkGeo = 'benchmarks/open/open.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end