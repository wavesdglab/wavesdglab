function [solA, matA, rhsA] = computeSolNum2D_UDG1(mesh, dofm, tau)

fprintf('Solver  : Call computeSolNumUDG1\n');

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

dofTRI = dofm.numDofTRI;
dofLIN = dofm.numDofLIN;

matA = [
    -1i*k*matM               -matDX                   -matDY                   sparse(dofTRI,2*dofLIN)   ;
    -matDX                   -1i*k*matM               sparse(dofTRI,dofTRI)    sparse(dofTRI,2*dofLIN)   ;
    -matDY                   sparse(dofTRI,dofTRI)    -1i*k*matM               sparse(dofTRI,2*dofLIN)   ;
    sparse(2*dofLIN,dofTRI)  sparse(2*dofLIN,dofTRI)  sparse(2*dofLIN,dofTRI)  sparse(1:2*dofLIN,1:2*dofLIN,1) ];

rhsA = [
    -1/(1i*k)*matM*solF ;
    zeros(dofTRI,1)     ;
    zeros(dofTRI,1)     ;
    zeros(2*dofLIN,1)   ];

% -------------------------------------------------------------------------
% Build surface terms
% -------------------------------------------------------------------------
fprintf('Solver  : Build surface terms\n');

dofLocTri      = [ 1 2 ; 2 3 ; 3 1 ];
dofLocTriNeigh = [ 2 1 ; 3 2 ; 1 3 ];

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
        idIntP = 0*dofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntU = 1*dofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntV = 2*dofTRI + dofm.locToGloTRI(tri,dofInt);
        
        % Global ID for (interior) edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        idIntS = 3*dofTRI + dofm.locToGloLIN(edgGlo,:);
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
        
        % Interior contributions
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau           * matMel;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5     * nx      * matMel;
        matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5     * ny      * matMel;
        matA(idIntP,idExtG) = matA(idIntP,idExtG) - 0.5               * matMel;
        
        matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5          * nx * matMel;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * nx * nx * matMel;
        matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5/tau * nx * ny * matMel;
        matA(idIntU,idExtG) = matA(idIntU,idExtG) + 0.5/tau      * nx * matMel;
        
        matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5          * ny * matMel;
        matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5/tau * nx * ny * matMel;
        matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5/tau * ny * ny * matMel;
        matA(idIntV,idExtG) = matA(idIntV,idExtG) + 0.5/tau      * ny * matMel;
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = 0*dofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*dofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*dofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            
            matA(idExtG,idExtP) = matA(idExtG,idExtP) - tau * sparse(1:2,1:2,1);
            matA(idExtG,idExtU) = matA(idExtG,idExtU) + nx  * sparse(1:2,1:2,1);
            matA(idExtG,idExtV) = matA(idExtG,idExtV) + ny  * sparse(1:2,1:2,1);
            
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
                    matA(idExtG,idIntP) = matA(idExtG,idIntP) + tau * sparse(1:2,1:2,1);
                    matA(idExtG,idIntU) = matA(idExtG,idIntU) + nx  * sparse(1:2,1:2,1);
                    matA(idExtG,idIntV) = matA(idExtG,idIntV) + ny  * sparse(1:2,1:2,1);
                    rhsA(idExtG) = rhsA(idExtG) + 2*tau*sol(idIntP);
                case 'NEU'
                    matA(idExtG,idIntP) = matA(idExtG,idIntP) - tau * sparse(1:2,1:2,1);
                    matA(idExtG,idIntU) = matA(idExtG,idIntU) - nx  * sparse(1:2,1:2,1);
                    matA(idExtG,idIntV) = matA(idExtG,idIntV) - ny  * sparse(1:2,1:2,1);
                    rhsA(idExtG) = rhsA(idExtG) - 2*(nx*solU(idIntP) + ny*solV(idIntP));
                case 'ABC'
                    
                    % Fix interior contributions
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*(1-tau)           * matMel;
                    
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5*(1-1/tau) * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matA(idIntU,idExtG) = matA(idIntU,idExtG) + 0.5*(1-1/tau)      * nx * matMel;
                    
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5*(1-1/tau) * ny * ny * matMel;
                    matA(idIntV,idExtG) = matA(idIntV,idExtG) + 0.5*(1-1/tau)      * ny * matMel;
                    
                    rhsA(idExtG) = rhsA(idExtG) + sol(idIntP) - (nx*solU(idIntP) + ny*solV(idIntP));
                otherwise
                    warning('Error - Bad BC.')
            end
        end
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

fprintf('Solver  : Solve ... \n');
solA = matA\rhsA;
solA = solA(1:dofTRI);

fprintf('---------------------------------------------------------\n');

end