function mesh = myBenchmark()

global h edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdisk Rdom Rpml

global PML_TYPE
if(isempty(PML_TYPE))
    PML_TYPE = 'Circular';
end

switch PML_TYPE
    case 'Rectangular'
        
        if(isempty(LdomX)) LdomX = 1.2; end
        if(isempty(LdomY)) LdomY = 1.2; end
        if(isempty(LpmlX)) LpmlX = 0.2; end
        if(isempty(LpmlY)) LpmlY = 0.2; end
        if(isempty(Rdisk)) Rdisk = 1; end
        
        % BCPML, BCObstacle
        edgTag = {201, 202};
        BC = {'NEU', 'NEU0'};
        edgTagToBC = containers.Map(edgTag,BC);
        
        linkMsh = 'output/mesh.msh';
        linkGeo = 'benchmarks/scattering_disk/scatteringDisk_RectangularPML.geo';
        system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ...
            ' -setnumber LdomX ' num2str(LdomX) ' -setnumber LdomY ' num2str(LdomY) ...
            ' -setnumber LpmlX ' num2str(LpmlX) ' -setnumber LpmlY ' num2str(LpmlY) ...
            ' -setnumber Rdisk ' num2str(Rdisk)]);
        mesh = readMesh2D(linkMsh);
        
    case 'Circular'
        
        if(isempty(Rdisk)) Rdisk = 1; end
        if(isempty(Rdom)) Rdom = 1.2; end
        if(isempty(Rpml)) Rpml = 0.2; end
        
        % BCPML, BCObstacle
        edgTag = {201, 202};
        BC = {'NEU', 'NEU0'};
        edgTagToBC = containers.Map(edgTag,BC);
        
        linkMsh = 'output/mesh.msh';
        linkGeo = 'benchmarks/scattering_disk/scatteringDisk_CircularPML.geo';
        system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ...
            ' -setnumber Rdisk ' num2str(Rdisk) ...
            ' -setnumber Rdom ' num2str(Rdom) ...
            ' -setnumber Rpml ' num2str(Rpml)]);
        mesh = readMesh2D(linkMsh);
        
end

end