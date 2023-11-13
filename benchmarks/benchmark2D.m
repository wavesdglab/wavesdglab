function mesh = benchmark2D(tag,h)

global BCWest BCNorth BCEast BCSouth TAGbench BCObstacle L L_PML R_disk BCPML l

TAGbench = tag;
switch tag
    case 'cavity'
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
    case 'scatteringPML'
        BCPML = 'DIR';
        BCObstacle = 'NEU';
        linkMsh = 'benchmarks/scattering/scatteringPML.msh';
        linkGeo = 'benchmarks/scattering/scatteringPML.geo';
    case 'scattering_rect'
        BCPML = 'DIR';
        BCObstacle = 'NEU';
        linkMsh = 'benchmarks/scattering_cav/rec_cavity.msh';
        linkGeo = 'benchmarks/scattering_cav/rec_cavity.geo';
    case 'scattering_square'
        BCPML = 'DIR';
        BCObstacle = 'NEU';
        linkMsh = 'benchmarks/scattering_cav/square_cavity.msh';
        linkGeo = 'benchmarks/scattering_cav/square_cavity.geo';
    otherwise
        warning('Error - No valid benchmark has been set.')
end

% disp(['L = ' L ', L_PML = ' L_PML ', R_disk = ' R_disk]);

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber L ' num2str(L) ' -setnumber L_PML ' num2str(L_PML) ' -setnumber R_disk ' num2str(R_disk) ' -setnumber l ' num2str(l)]);
mesh = readMesh2D(linkMsh);

end