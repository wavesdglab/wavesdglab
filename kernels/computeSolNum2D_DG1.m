function [solA, matA, rhsA] = computeSolNum2D_DG1(mesh, dofm, tau)

fprintf('Solver  : Call computeSolNumDG1\n');

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

matA = [
    -1i*k*matM  -matDX                       -matDY                      ;
    -matDX      -1i*k*matM                   sparse(numDofTRI,numDofTRI) ;
    -matDY      sparse(numDofTRI,numDofTRI)  -1i*k*matM                  ];

rhsA = [
    -1/(1i*k)*matM*solF ;
    zeros(numDofTRI,1) ;
    zeros(numDofTRI,1) ];

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
        
        % Elemental matrices
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        verEdg = mesh.mapEdgToVer(edgGlo,:);
        V1 = mesh.coord(verEdg(1),:);
        V2 = mesh.coord(verEdg(2),:);
        [matMel, ~, ~] = buildMatrixElemLIN(V1,V2,dofm.degree);
        
        % Global ID for interior unknowns
        dofInt = dofLocTri(fac,:);
        idIntP = dofm.locToGloTRI(tri,dofInt);
        idIntU = idIntP + numDofTRI;
        idIntV = idIntU + numDofTRI;
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = idExtP + numDofTRI;
            idExtV = idExtU + numDofTRI;
            
            matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau         * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5*nx          * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5*ny          * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau         * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5*nx          * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + 0.5*ny          * matMel;
            
            matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5*nx          * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * nx*nx * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5/tau * nx*ny * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5*nx          * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * nx*nx * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) - 0.5/tau * nx*ny * matMel;
            
            matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5*ny          * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5/tau * nx*ny * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5/tau * ny*ny * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + 0.5*ny          * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) - 0.5/tau * nx*ny * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) - 0.5/tau * ny*ny * matMel;
            
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
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau             * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx              * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny              * matMel;
                    
                    rhsA(idIntP) = rhsA(idIntP) + matMel * tau*sol(idIntP);
                    rhsA(idIntU) = rhsA(idIntU) - matMel * nx*sol(idIntP);
                    rhsA(idIntV) = rhsA(idIntV) - matMel * ny*sol(idIntP);
                    
                case 'NEU'
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + nx              * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/tau * nx*nx   * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 1/tau * nx*ny   * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + ny              * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 1/tau * nx*ny   * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 1/tau * ny*ny   * matMel;
                    
                    rhsA(idIntP) = rhsA(idIntP) - matMel * (nx*solU(idIntP) + ny*solV(idIntP));
                    rhsA(idIntU) = rhsA(idIntU) + matMel * nx*solU(idIntP);
                    rhsA(idIntV) = rhsA(idIntV) + matMel * ny*solV(idIntP);
                    
                case 'ABC'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5             * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5*nx          * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5*ny          * matMel;
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5*nx          * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5     * nx*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5     * nx*ny * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5*ny          * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5     * nx*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5     * ny*ny * matMel;
                    
                    gchar = sol(idIntP) - (nx*solU(idIntP) + ny*solV(idIntP));
                    rhsA(idIntP) = rhsA(idIntP) + matMel * 0.5 * gchar   ;
                    rhsA(idIntU) = rhsA(idIntU) - matMel * 0.5 * nx*gchar;
                    rhsA(idIntV) = rhsA(idIntV) - matMel * 0.5 * ny*gchar;
                    
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
solA = solA(1:numDofTRI);

fprintf('---------------------------------------------------------\n');

end