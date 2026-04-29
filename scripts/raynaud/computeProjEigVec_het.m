% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the nbEigVec closest eigenvectors to k for the Laplacian problem on a rectagular domain with homogeneous Neumann boundary conditions

% Let u (resp. u_PML) be the analytical expression of a quasimode without (resp. with) a complex stretching in the PML region, M the mass matrix without PML and M_PML the mass matrix with PML.
% eigenvec : the projection of u on the finite element space (without M_PML)
% eigenvec_PML : the projection of u_PML on the finite element space with M_PML
% y = M*eigenvec
% y_PML = M_PML*eigenvec_PML

function [nbEigVec,eigenvec,eigenvec_PML,y,y_PML] = computeProjEigVec_het(mesh, dofm, nbEigVec)

global edgTagToBC Options cObj Rdom Rpml k

noptique = 1/cObj;

if nbEigVec == 0
    eigenvec = [];
    eigenvec_PML = [];
    y = [];
    y_PML = [];
    return;
end

mn = computeCloseEigVec_het(nbEigVec);

nbEigVec = 2*size(mn, 1);

% Quadrature and shape functions
degreeQ = 2*dofm.degree;
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
shapeQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);

% Build matrix and RHS vector
matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
matP_PML = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, 2*size(mn, 1));
rhsP_PML = zeros(dofm.numDofTRI, 2*size(mn, 1));
for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % Reference solution
    refQ = zeros(size(xQ,1), 2*size(mn, 1));
    refQ_PML = zeros(size(xQ,1), 2*size(mn, 1));

    for i=1:size(mn, 1)
        m = mn(i, 1);
        n = mn(i, 2);
        aj = -airy_zero(n+1);
        % si la moyenne des coordoonées du triangle est <= 1 : refQ =  exp(1i * m * theta(x,y)) * (gslsfairyAi(-aj - 2.^(1./3.)*sig(h), 0.) + h*2.^(1./3.)*noptique*gslsfairyAideriv(-aj - 2.^(1./3.)*sig(h), 0.)/sqrt(noptique^2-1.) )
        % sinon : refQ = exp(1i * m * theta(x,y)) * 2.^(1./3.) * h* gslsfairyAideriv(-aj,0) * noptique / (sqrt(noptique^2-1.)) * exp(- rho(h)/noptique * sqrt(noptique^2-1.))
        meanX = mean([V1(1), V2(1), V3(1)]);
        meanY = mean([V1(2), V2(2), V3(2)]);
        if sqrt(meanX^2 + meanY^2) <= 1
            refQ(:,2*i-1) = exp(1i * m * atan2(yQ, xQ)) .* (airy(-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) + (1./m).^(1/3)*2.^(1./3.)*noptique*airy(1,-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) ./ sqrt(noptique^2-1.) );
            refQ(:,2*i) = exp(-1i * m * atan2(yQ, xQ)) .* (airy(-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) + (1./m).^(1/3)*2.^(1./3.)*noptique*airy(1,-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) ./ sqrt(noptique^2-1.) );
            refQ_PML(:,2*i-1) = exp(1i * m * atan2(yQ, xQ)) .* (airy(-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) + (1./m).^(1/3)*2.^(1./3.)*noptique*airy(1,-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) ./ sqrt(noptique^2-1.) );
            refQ_PML(:,2*i) = exp(-1i * m * atan2(yQ, xQ)) .* (airy(-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) + (1./m).^(1/3)*2.^(1./3.)*noptique*airy(1,-aj - 2.^(1./3.)*sigma(xQ,yQ,m)) ./ sqrt(noptique^2-1.) );
        else
            refQ(:,2*i-1) = exp(1i * m * atan2(yQ, xQ)) .* 2.^(1./3.) * (1./m).^(1/3) * airy(1,-aj) * noptique / (sqrt(noptique^2-1.)) .* exp(- rho(xQ,yQ,m)/noptique * sqrt(noptique^2-1.));
            refQ(:,2*i) = exp(-1i * m * atan2(yQ, xQ)) .* 2.^(1./3.) * (1./m).^(1/3) * airy(1,-aj) * noptique / (sqrt(noptique^2-1.)) .* exp(- rho(xQ,yQ,m)/noptique * sqrt(noptique^2-1.));
            refQ_PML(:,2*i-1) = exp(1i * m * atan2(yQ, xQ)) .* 2.^(1./3.) * (1./m).^(1/3) * airy(1,-aj) * noptique / (sqrt(noptique^2-1.)) .* exp(- rhoPML(xQ,yQ,m)/noptique * sqrt(noptique^2-1.));
            refQ_PML(:,2*i) =  exp(-1i * m * atan2(yQ, xQ)) .* 2.^(1./3.) * (1./m).^(1/3) * airy(1,-aj) * noptique / (sqrt(noptique^2-1.)) .* exp(- rhoPML(xQ,yQ,m)/noptique * sqrt(noptique^2-1.));
        end
        
    end
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if ~strcmp(Options.Basis,'Lagrange')
        if(ver(1) > ver(2))
            orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(2) > ver(3))
            orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(3) > ver(1))
            orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);

    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    
    % Local matrix and RHS vector
    weightsQ = weightsTriQ .* detJdxdu;
    matPel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
    matP_PMLel = matPel;
    rhsPel = transpose(shapeOrQ) * (weightsQ .* refQ);
    rhsP_PML = transpose(shapeOrQ) * (weightsQ .* refQ_PML);
    if(~isempty(Rdom))
        rQ = sqrt(xQ.*xQ + yQ.*yQ);
        if (mean(rQ) >= Rdom)
            % cosT = xQ./rQ;
            % sinT = yQ./rQ;
            sigmaPml = 1./(Rpml-(rQ-Rdom));
            sigmaPmlInt = -log(1-(rQ-Rdom)/Rpml);
            gammaPmlR = ones(size(rQ)) - sigmaPml/(1i*k(tri));
            gammaPmlT = ones(size(rQ)) - sigmaPmlInt/(1i*k(tri))./rQ;
            % invJacXX = (1./gammaPmlR) .* cosT.*cosT + (1./gammaPmlT) .* (sinT.*sinT);
            % invJacXY = (1./gammaPmlR) .* cosT.*sinT - (1./gammaPmlT) .* (cosT.*sinT);
            % invJacYY = (1./gammaPmlR) .* sinT.*sinT + (1./gammaPmlT) .* (cosT.*cosT);
            detJdxdu = detJdxdu .* gammaPmlR .* gammaPmlT;
            % shapeDxQnew = invJacXX .* shapeDxQ + invJacXY .* shapeDyQ;
            % shapeDyQnew = invJacXY .* shapeDxQ + invJacYY .* shapeDyQ;
            % shapeDxQ = shapeDxQnew;
            % shapeDyQ = shapeDyQnew;

            % Elemental matrices in PML
            weightsQ = weightsTriQ .* detJdxdu;
            matMpml = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
            matP_PMLel = matMpml;
        end
    end
    
    % Assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    matP_PML(dof,dof) = matP_PML(dof,dof) + matP_PMLel;
    rhsP(dof,:) = rhsP(dof,:) + rhsPel;
    rhsP_PML(dof,:) = rhsP_PML(dof,:) + rhsP_PML;
end


dofDIR = [];
cacheDIR = zeros(dofm.numDofTRI);
for edgBnd=1:mesh.numEdgBnd
    dof = dofm.locToGloBND(edgBnd,:);
    % Boundary condition
    switch edgTagToBC(mesh.tagEdgBnd(edgBnd))
        case 'DIR0'
            dofDIR = [dofDIR ; dof];
            cacheDIR(dof) = zeros(size(dof,1),size(dof,2));
    end
end
if(~isempty(dofDIR))
    dofDIR = unique(dofDIR);
    matP(:,dofDIR) = 0;
    matP(dofDIR,:) = 0;
    matP(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
    rhsP(dofDIR,:) = 0;
    rhsP_PML(dofDIR,:) = 0;
    matP_PML(:,dofDIR) = 0;
    matP_PML(dofDIR,:) = 0;
    matP_PML(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end


% Solution
eigenvec = matP\rhsP;
eigenvec_PML = matP_PML\rhsP_PML;
y = rhsP;
y_PML = rhsP_PML;


% for i=1:size(mn, 1)
%     m = mn(i, 1);
%     n = mn(i, 2);
%     filename = 'output/eigenvec' + string(m)+'*'+string(n) + '.pos';
%     fieldname = 'eigenvec' + string(m)+'*'+string(n);
%     writeField2D(dofm, mesh, eigenvec(:,i), filename, fieldname);
%     fprintf('field [%i] saved \n', i);
%     disp(['m = ' num2str(m) ', n = ' num2str(n) ', k_jm = ' num2str(kjm(m, n))]);
% end

end

function sig = sigma(x,y,m)
    h = (1./m).^(1/3);
    sig = h.^(-2).*(sqrt(x.^2 + y.^2) - 1);
end

function val = rho(x,y,m)
    h = (1./m).^(1/3);
    val = h.^(-3).*(sqrt(x.^2 + y.^2) - 1);
end

function val = rhoPML(x,y,m)
    global omega Rdom Rpml
    h = (1./m).^(1/3);
    r = sqrt(x.^2 + y.^2);
    if mean(r) >= Rdom
        rprime = r - 1/(1i*omega) * log((Rpml - Rdom) ./(max(1e-10, Rpml - r)));
    else
        rprime = r;
    end
    val = h.^(-3).*(rprime - 1);
end

function z = airy_zero(k)
    Ai = @(x) airy(0,x);
    x0 = -((3*pi/2)*(k - 1/4)).^(2/3);
    z = zeros(size(k));
    for i = 1:numel(k)
        z(i) = fzero(Ai, [x0(i) - 1, x0(i) + 1]);
    end
end

function karray = coeffsk(j)
    global cObj
    noptique = 1/cObj;
    
    aj = -airy_zero(j+1);

    karray = zeros(10, numel(j));
    karray(1, :) = aj / 2.0;
    karray(2, :) = -noptique / (2 * sqrt(noptique^2 - 1));
    karray(3, :) = 3 * aj.^2 / 40.0;
    karray(4, :) = -aj .* noptique * (3*noptique^2 - 2*noptique^2) ...
           / (12 * (noptique^2 - 1)^(3/2));
    karray(5, :) = 0.125 * (1/35 - aj.^3 / 350);
    karray(6, :) = - (aj.^2 .* noptique * (-noptique^4 + 4*noptique^2)) ...
           / (80 * (noptique^2 - 1)^(5/2));
    karray(7, :) = -aj / 144.0 .* (1/175 + 479*aj.^3 / 7000 ...
           - noptique * (-2*noptique^5) / (noptique^2 - 1)^3);
    karray(8, :) = -noptique^3 / (35 * (noptique^2 - 1)^(7/2)) .* ...
           (4*aj.^3/5 + 13/64 ...
           - (10 - aj.^3) * noptique^2 * (noptique^2 - 4) / 160);
    karray(9, :) = aj.^2 / 240 .* ...
           (20231*aj.^3 / 1078000 - 551/10780 ...
           + noptique^6 * (noptique^2 - 7) / (noptique^2 - 1)^4);
    karray(10, :) = -aj .* noptique^3 / (1260 * (noptique^2 - 1)^(9/2)) .* ...
            ( noptique^6 * (479*aj.^3/960 + 353/72) ...
            - noptique^4 * (479*aj.^3/200 + 0.2) ...
            + noptique^2 * (1009*aj.^3/100 + 4919/160) ...
            + 821*aj.^3/75 + 1847/240 );
end

function val = kjm(m, j)
    global cObj
    noptique = 1/cObj;
    kcoeffs = coeffsk(j(:,1));
    i = 0:length(kcoeffs)-1;
    sum_val = ones(size(m));
    for idx = 1:size(m, 1)
        for jdx = 1:size(m, 2)
            for kdx = 1:length(i)
                sum_val(idx,jdx) = sum_val(idx,jdx) + (2/m(idx,jdx))^(2/3) * kcoeffs(kdx, idx) * (2/m(idx,jdx))^(i(kdx)/3);
            end
        end
    end
    val = (m ./ noptique) .* sum_val;
end


function indices = computeCloseEigVec_het(nb)
    global omega
    limit = 100;
    M = 1:limit;
    N = 0:2;
    [M, N] = meshgrid(M, N);

    quasi_resonances = kjm(M, N);

    diff = abs(quasi_resonances - omega);

    [~, sorted_indices] = sort(diff(:));


    indices = sorted_indices(1:nb);


    [n, m] = ind2sub(size(diff), indices);

    n = n-1;

    indices = [m(1:nb), n(1:nb)];
end
