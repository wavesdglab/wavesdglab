function [solA, matA, rhsA] = computeSolNum2D_UDG2(mesh, dofm, tau)

fprintf('Solver  : Call computeSolNumUDG2\n');

global k BCWest BCNorth BCEast BCSouth

% -------------------------------------------------------------------------
% Get data
% -------------------------------------------------------------------------

fprintf('Solver  : Get data\n');
x = mesh.coord(:,1);
y = mesh.coord(:,2);
x = x(mesh.mapTriToVer)';
y = y(mesh.mapTriToVer)';
x = x(:);
y = y(:);
[sol, solDx, solDy, solF] = mySol(x,y);

% -------------------------------------------------------------------------
% Build volume terms
% -------------------------------------------------------------------------
fprintf('Solver  : Build volume terms\n');

[matM, matK, ~, ~] = buildMatrixGlo2D_DG(mesh, dofm);

dofTRI = dofm.numDofTRI;
dofLIN = dofm.numDofLIN;

matA = [
    matK - k^2*matM          sparse(dofTRI,2*dofLIN)         ;
    sparse(2*dofLIN,dofTRI)  sparse(1:2*dofLIN,1:2*dofLIN,1) ];

rhsA = [
    matM * solF       ;
    zeros(2*dofLIN,1) ];

% -------------------------------------------------------------------------
% Build surface terms
% -------------------------------------------------------------------------
fprintf('Solver  : Build surface terms\n');

dofLocTri      = [ 1 2 ; 2 3 ; 3 1 ];
dofLocTriNeigh = [ 2 1 ; 3 2 ; 1 3 ];

dofDIR = [];

for tri=1:mesh.numTri
    
    % Compute exterior normals
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    normal = getNormalTRI(V1,V2,V3);
    
    % Loop over faces
    for fac = 1:3
        
        % Global ID for interior unknowns
        dofInt = dofLocTri(fac,:);
        idIntP = dofm.locToGloTRI(tri,dofInt);
        
        % Global ID for (interior) edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        idIntS = dofTRI + dofm.locToGloLIN(edgGlo,:);
        if(mesh.mapTriToEdg(tri,fac) > 0)
            idIntG = idIntS;
            idExtG = idIntS + dofLIN;
        else
            idIntG = idIntS([2 1]) + dofLIN;
            idExtG = idIntS([2 1]);
        end
        
        % Elemental matrices
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        verEdg = mesh.listEdg(edgGlo,:);
        V1 = mesh.coord(verEdg(1),:);
        V2 = mesh.coord(verEdg(2),:);
        [matMel, ~, ~] = buildMatrixElemLIN(V1,V2,dofm.degree);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            matA(idIntP,idIntP) = matA(idIntP,idIntP) - 1i*tau*k * matMel;
            matA(idIntP,idIntG) = matA(idIntP,idIntG) -            matMel;
            
            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            
            matA(idIntG,idExtG) = matA(idIntG,idExtG) +            matMel;
            matA(idIntG,idIntP) = matA(idIntG,idIntP) + 2i*tau*k * matMel;
            
        else
            
            switch mesh.tagEdg(edgGlo)
                case 1
                    BC = BCWest;
                case 2
                    BC = BCNorth;
                case 3
                    BC = BCEast;
                case 4
                    BC = BCSouth;
                otherwise
                    warning('Error - Bad TAG.')
            end
            
            switch BC
                case 'DIR'
                    dofDIR = [dofDIR ; idIntP];
                case 'NEU'
                    rhsA(idIntP) = rhsA(idIntP) + matMel * (nx*solDx(idIntP) + ny*solDy(idIntP));
                case 'ABC'
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) - 1i*k * matMel;
                    rhsA(idIntP) = rhsA(idIntP) + matMel * (nx*solDx(idIntP) + ny*solDy(idIntP) - 1i*k*sol(idIntP));
                otherwise
                    warning('Error - Bad BC.')
            end
        end
    end
end

if(~isempty(dofDIR))
    dofDIR = unique(dofDIR);
    rhsA = rhsA - matA(:,dofDIR)*sol(dofDIR);
    rhsA(dofDIR) = sol(dofDIR);
    matA(dofDIR,:) = 0;
    matA(:,dofDIR) = 0;
    matA(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

fprintf('Solver  : Solve ... \n');
solA = matA\rhsA;
solA = solA(1:dofTRI);

fprintf('---------------------------------------------------------\n');

end