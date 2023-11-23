function mesh = benchmark2D(tag,h)

global BCWest BCNorth BCEast BCSouth TAGbench

TAGbench = tag;
switch tag
    case 'cavity'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/cavity/cavity.msh';
        linkGeo = 'benchmarks/cavity/cavity.geo';
    case 'cavity (heterogeneous)'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/cavity_heterogeneous/cavity_heterogeneous.msh';
        linkGeo = 'benchmarks/cavity_heterogeneous/cavity_heterogeneous.geo';
    case 'waveguide'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'ABC';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/waveguide/waveguide.msh';
        linkGeo = 'benchmarks/waveguide/waveguide.geo';
    case 'waveguide (heterogeneous)'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'ABC';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/waveguide_heterogeneous/waveguide_heterogeneous.msh';
        linkGeo = 'benchmarks/waveguide_heterogeneous/waveguide_heterogeneous.geo';    
    case 'open'
        BCWest  = 'ABC';
        BCNorth = 'ABC';
        BCEast  = 'ABC';
        BCSouth = 'ABC';
        linkMsh = 'benchmarks/open/open.msh';
        linkGeo = 'benchmarks/open/open.geo';
    case 'open (heterogeneous)'
        BCWest  = 'ABC';
        BCNorth = 'ABC';
        BCEast  = 'ABC';
        BCSouth = 'ABC';
        linkMsh = 'benchmarks/open_heterogeneous/open_heterogeneous.msh';
        linkGeo = 'benchmarks/open_heterogeneous/open_heterogeneous.geo';
    case 'square'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/square/square.msh';
        linkGeo = 'benchmarks/square/square.geo';
    otherwise
        warning('Error - No valid benchmark has been set.')
end

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end