function buildMatrixGlo2D_HDGsave()

global Coord
global NumEdg
global NumTri
global ListTri ListEdg
global mapTriToEdg

numDofSur = 2*NumEdg;
numDofVol = 3*NumTri;

% -------------------------------------------------------------------------
% Compute surface-to-surface matrices
% -------------------------------------------------------------------------

global IntMM
IntMM = sparse(numDofSur,numDofSur);   % Mass matrix

for edg=1:NumEdg
    
    % Elemental matrices
    [Mel] = buildMatrixElemLIN(edg);
    
    % Matrix assembling
    dofGloEdg = 2*(edg-1)+1:2*edg;
    IntMM(dofGloEdg,dofGloEdg) = Mel;
end


% -------------------------------------------------------------------------
% Compute volume-to-volume matrices + volume-to-surface matrices
% -------------------------------------------------------------------------

global VolMM VolKK VolDX VolDY VolNX VolNY VolSS

vecI  = zeros(9*NumTri,1);
vecJ  = zeros(9*NumTri,1);
vecMM = zeros(9*NumTri,1);   % Mass matrix
vecKK = zeros(9*NumTri,1);   % Stiffness matrix
vecDX = zeros(9*NumTri,1);   % Differentiation matrix (x)
vecDY = zeros(9*NumTri,1);   % Differentiation matrix (y)
vecNX = zeros(9*NumTri,1);   % Element-normalX mass matrix
vecNY = zeros(9*NumTri,1);   % Element-normalY mass matrix
vecSS = zeros(9*NumTri,1);   % Element-surface mass matrix

global NX2I NY2I SS2I
NX2I = sparse(numDofSur,numDofVol);
NY2I = sparse(numDofSur,numDofVol);
SS2I = sparse(numDofSur,numDofVol);

vecIel = [1 1 1 ; 2 2 2 ; 3 3 3 ];
vecJel = [1 2 3 ; 1 2 3 ; 1 2 3 ];

for tri=1:NumTri
    
    % Elemental matrices
    [Mel, Kel, DXel, DYel] = buildMatrixElemTRI(tri);
    Mel = [1 2 3; 4 5 6; 7 8 9];
    Kel = [1 2 3; 4 5 6; 7 8 9];
    DXel = [1 2 3; 4 5 6; 7 8 9];
    DYel = [1 2 3; 4 5 6; 7 8 9];
    
    % Matrix assembling
    index = 9*(tri-1)+1:9*tri;
    vecI(index) = 3*(tri-1)+vecIel(:);
    vecJ(index) = 3*(tri-1)+vecJel(:);
    vecMM(index) = Mel(:);
    vecKK(index) = Kel(:);
    vecDX(index) = DXel(:);
    vecDY(index) = DYel(:);
    
    % Compute exterior normals
    m = ListTri(tri,:);
    V1 = Coord(m(1),:);
    V2 = Coord(m(2),:);
    V3 = Coord(m(3),:);
    N1 = -(V3-V1) + ((V3-V1)*(V2-V1)') * (V2-V1)/norm(V2-V1)^2; N1 = N1/norm(N1);
    N2 = -(V1-V2) + ((V1-V2)*(V3-V2)') * (V3-V2)/norm(V3-V2)^2; N2 = N2/norm(N2);
    N3 = -(V2-V3) + ((V2-V3)*(V1-V3)') * (V1-V3)/norm(V1-V3)^2; N3 = N3/norm(N3);
    N = [N1 ; N2 ; N3];
    
    %     hold off
    %     plot([V1(1) V2(1) V3(1) V1(1)], [V1(2) V2(2) V3(2) V1(2)]);
    %     hold on
    %     plot([(V1(1)+V2(1))/2 (V1(1)+V2(1))/2+N1(1)], [(V1(2)+V2(2))/2 (V1(2)+V2(2))/2+N1(2)]);
    %     plot([(V2(1)+V3(1))/2 (V2(1)+V3(1))/2+N2(1)], [(V2(2)+V3(2))/2 (V2(2)+V3(2))/2+N2(2)]);
    %     plot([(V3(1)+V1(1))/2 (V3(1)+V1(1))/2+N3(1)], [(V3(2)+V1(2))/2 (V3(2)+V1(2))/2+N3(2)]);
    %     axis([-1 1 -1 1]);
    %     pause;
    
    % Indices
    edgGlo = mapTriToEdg(tri,:);
    dofGloTri = 3*(tri-1)+1:3*tri;
    dofGloEdg = [ 2*(abs(edgGlo(1))-1)+1 2*abs(edgGlo(1)) ;
                  2*(abs(edgGlo(2))-1)+1 2*abs(edgGlo(2)) ;
                  2*(abs(edgGlo(3))-1)+1 2*abs(edgGlo(3)) ];
    
    nodLocFacToNodLocTri = [ 1 2 ; 2 3 ; 3 1 ]; % nodLocEdg to nodLocTri
    for f = 1:3
        if (edgGlo(f) < 0)
            tmp = nodLocFacToNodLocTri(f,:);
            nodLocFacToNodLocTri(f,:) = [tmp(2) tmp(1)];
        end
        if ((m(nodLocFacToNodLocTri(f,1)) ~= ListEdg(abs(edgGlo(f)),1)) || ...
            (m(nodLocFacToNodLocTri(f,2)) ~= ListEdg(abs(edgGlo(f)),2)))
             fprintf("(%i %i) (%i %i)\n", m(nodLocFacToNodLocTri(f,1)), m(nodLocFacToNodLocTri(f,2)), ...
                      ListEdg(abs(edgGlo(f)),1), ListEdg(abs(edgGlo(f)),2));
        end
    end
    
    
    % Compute elemental matrices
    NXel = zeros(3,3);
    NYel = zeros(3,3);
    SSel = zeros(3,3);
    for f = 1:3
        [Mel, ~, ~] = buildMatrixElemLIN(abs(edgGlo(f)));
        nodLocTri = nodLocFacToNodLocTri(f,:);
        
        NXel(nodLocTri, nodLocTri) = NXel(nodLocTri, nodLocTri) + N(f,1)*Mel;
        NYel(nodLocTri, nodLocTri) = NYel(nodLocTri, nodLocTri) + N(f,2)*Mel;
        SSel(nodLocTri, nodLocTri) = SSel(nodLocTri, nodLocTri) + Mel;
        
        NX2I(dofGloEdg(f,:), dofGloTri(nodLocTri)) = N(f,1)*Mel;
        NY2I(dofGloEdg(f,:), dofGloTri(nodLocTri)) = N(f,2)*Mel;
        SS2I(dofGloEdg(f,:), dofGloTri(nodLocTri)) = Mel;
    end
    
    % Matrix assembling
    vecNX(index) = NXel(:);
    vecNY(index) = NYel(:);
    vecSS(index) = SSel(:);
end

VolMM = sparse(vecI,vecJ,vecMM);   % Mass matrix
VolKK = sparse(vecI,vecJ,vecKK);   % Stiffness matrix
VolDX = sparse(vecI,vecJ,vecDX);   % Differentiation matrix (x)
VolDY = sparse(vecI,vecJ,vecDY);   % Differentiation matrix (y)
VolNX = sparse(vecI,vecJ,vecNX);   % Element-normalX mass matrix
VolNY = sparse(vecI,vecJ,vecNY);   % Element-normalY mass matrix
VolSS = sparse(vecI,vecJ,vecSS);   % Element-surface mass matrix

end