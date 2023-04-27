function mesh = benchmark2D(tag,h)

global BCWest BCNorth BCEast BCSouth TAGbench BCCircle L L_PML R_disk

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
    case 'scatteringPML'
        BCWest = 'DIR';
        BCNorth = 'DIR';
        BCEast = 'DIR';
        BCSouth = 'DIR';
        BCCircle = 'NEU';
        linkMsh = 'benchmarks/scattering/scatteringPML.msh';
        linkGeo = 'benchmarks/scattering/scatteringPML.geo';

    otherwise
        warning('Error - No valid benchmark has been set.')
end

% disp(['L = ' L ', L_PML = ' L_PML ', R_disk = ' R_disk]);

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber L ' num2str(L) ' -setnumber L_PML ' num2str(L_PML) ' -setnumber R_disk ' num2str(R_disk)]);
mesh = readMesh2D(linkMsh);

end