function mesh = benchmark2D(tag,h)

global BCWest BCNorth BCEast BCSouth TAGbench BCCircle

TAGbench = tag;
switch tag
    case 'cavity'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/cavity/cavity.msh';
        linkGeo = 'benchmarks/cavity/cavity.geo';
    case 'cavity2'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'DIR';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/cavity/cavity.msh';
        linkGeo = 'benchmarks/cavity/cavity.geo';
    case 'waveguide'
        BCWest  = 'DIR';
        BCNorth = 'DIR';
        BCEast  = 'ABC';
        BCSouth = 'DIR';
        linkMsh = 'benchmarks/waveguide/waveguide.msh';
        linkGeo = 'benchmarks/waveguide/waveguide.geo';
    case 'open'
        BCWest  = 'ABC';
        BCNorth = 'ABC';
        BCEast  = 'ABC';
        BCSouth = 'ABC';
        linkMsh = 'benchmarks/open/open.msh';
        linkGeo = 'benchmarks/open/open.geo';
    case 'scatteringHard'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'NEU';
        linkMsh = 'benchmarks/scattering/scattering.msh';
        linkGeo = 'benchmarks/scattering/scattering.geo';
    case 'scatteringSoft'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'DIR';
        linkMsh = 'benchmarks/scattering/scattering.msh';
        linkGeo = 'benchmarks/scattering/scattering.geo';
    case 'scatteringPML_01'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'NEU';
        linkMsh = 'benchmarks/scattering/scatteringPML_01.msh';
        linkGeo = 'benchmarks/scattering/scatteringPML_01.geo';
    case 'scatteringPML_02'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'NEU';
        linkMsh = 'benchmarks/scattering/scatteringPML_02.msh';
        linkGeo = 'benchmarks/scattering/scatteringPML_02.geo';
    case 'scatteringPML_05'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'NEU';
        linkMsh = 'benchmarks/scattering/scatteringPML_05.msh';
        linkGeo = 'benchmarks/scattering/scatteringPML_05.geo';

    otherwise
        warning('Error - No valid benchmark has been set.')
end

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end