function [solA, matA, rhsA] = computeSolNum2D_HDG1(mesh, dofm, tau)

fprintf('Solver  : Call computeSolNumHDG1\n');

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
solU = solDx/(1i*k);
solV = solDy/(1i*k);

% -------------------------------------------------------------------------
% Build volume terms
% -------------------------------------------------------------------------
fprintf('Solver  : Build volume terms\n');

[matM, ~, matDX, matDY] = buildMatrixGlo2D_DG(mesh, dofm);

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;

matA = [
    -1i*k*matM                   -matDX                       -matDY                       sparse(numDofTRI,numDofLIN) ;
    -matDX                       -1i*k*matM                   sparse(numDofTRI,numDofTRI)  sparse(numDofTRI,numDofLIN) ;
    -matDY                       sparse(numDofTRI,numDofTRI)  -1i*k*matM                   sparse(numDofTRI,numDofLIN) ;
    sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofLIN) ];

rhsA = [
    -1/(1i*k)*matM*solF ;
    zeros(numDofTRI,1)  ;
    zeros(numDofTRI,1)  ;
    zeros(numDofLIN,1)  ];

% -------------------------------------------------------------------------
% Build surface terms
% -------------------------------------------------------------------------
fprintf('Solver  : Build surface terms\n');

dofLocTri = [ 1 2 ; 2 3 ; 3 1 ];

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
        idIntP = 0*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntU = 1*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntV = 2*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        idIntS = 3*numDofTRI + dofm.locToGloLIN(edgGlo,:);
        if(mesh.mapTriToEdg(tri,fac) < 0)
            tmp = idIntS;
            idIntS(1) = tmp(2);
            idIntS(2) = tmp(1);
        end
        
        % Elemental matrices
        verEdg = mesh.mapEdgToVer(edgGlo,:);
        V1 = mesh.coord(verEdg(1),:);
        V2 = mesh.coord(verEdg(2),:);
        [matMel, ~, ~] = buildMatrixElemLIN(V1,V2,dofm.degree);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Surface terms for the volume fields
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau * matMel;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx  * matMel;
        matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny  * matMel;
        matA(idIntP,idIntS) = matA(idIntP,idIntS) - tau * matMel;
        matA(idIntU,idIntS) = matA(idIntU,idIntS) + nx  * matMel;
        matA(idIntV,idIntS) = matA(idIntV,idIntS) + ny  * matMel;
        
        % Surface terms for the surface field
        triNeigh = mesh.mapTriToTri(tri,fac);
        if (triNeigh > 0)
            matA(idIntS,idIntP) = matA(idIntS,idIntP) + tau * matMel;
            matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
            matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
            matA(idIntS,idIntS) = matA(idIntS,idIntS) - tau * matMel;
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
                    warning('Error - Bad TAG.');
            end
            switch BC
                case 'DIR'
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) + matMel;
                    rhsA(idIntS) = rhsA(idIntS) + matMel * sol(idIntP);
                case 'NEU'
                    matA(idIntS,idIntP) = matA(idIntS,idIntP) + tau * matMel;
                    matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
                    matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) - tau * matMel;
                    rhsA(idIntS) = rhsA(idIntS) + matMel * (nx*solU(idIntP) + ny*solV(idIntP));
                case 'ABC'
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + (1-tau) * matMel;
                    matA(idIntP,idIntS) = matA(idIntP,idIntS) - (1-tau) * matMel;
                    % 
                    matA(idIntS,idIntP) = matA(idIntS,idIntP) +       matMel;
                    matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
                    matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) - 2   * matMel;
                    rhsA(idIntS) = rhsA(idIntS) + matMel * (nx*solU(idIntP) - ny*solV(idIntP) - sol(idIntP));
                otherwise
                    warning('Error - Bad BC.');
            end
        end
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

fprintf('Solver  : Solve ... \n');
solA = matA\rhsA;
solA = solA(1:numDofTRI);

fprintf('---------------------------------------------------------\n');

end