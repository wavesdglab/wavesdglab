function mesh = benchmark2D(tag,h)

global BCWest BCNorth BCEast BCSouth TAGbench

TAGbench = tag;
switch tag
    case 'cavity'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks2D/cavity/cavity.msh';
        linkGeo = 'benchmarks2D/cavity/cavity.geo';
    case 'waveguide'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'ABC';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks2D/waveguide/waveguide.msh';
        linkGeo = 'benchmarks2D/waveguide/waveguide.geo';
    case 'open'
        BCWest  = 'ABC';
        BCNorth = 'ABC';
        BCEast  = 'ABC';
        BCSouth = 'ABC';
        linkMsh = 'benchmarks2D/open/open.msh';
        linkGeo = 'benchmarks2D/open/open.geo';
    otherwise
        warning('Error - No valid benchmark has been set.')
end

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh(linkMsh);

end