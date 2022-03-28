% Build elemental matrices on the TRIANGLE

function [matM, matK, matDX, matDY] = buildMatrixElemTRI(V1, V2, V3, degree)

degreeQ = 2*degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
valQ = functionsShapeTRI(uQ, vQ, degree);
[valDx, valDy] = functionsShapeDerTRI(uQ, vQ, degree);

N = size(valQ,2);

% Reference line [-1,1]
matMref = zeros(N,N);
matKref = zeros(N,N);
matDXref = zeros(N,N);
matDYref = zeros(N,N);
for i=1:N
    for j=1:N
        matMref(i,j)  = weights' * (valQ(:,i) .* valQ(:,j));
        matKref(i,j)  = weights' * (valDx(:,i) .* valDx(:,j) + valDy(:,i) .* valDy(:,j));
        matDXref(i,j) = weights' * (valQ(:,i) .* valDx(:,j));
        matDYref(i,j) = weights' * (valQ(:,i) .* valDy(:,j));
    end
end

% Jacobian ([V1,V2,V3] <=> [(0,0),(0,-1),(0,1)])
J = [(V2-V1)' (V3-V1)'];
detJ = abs(det(J));
gradW = J'\[-1 1 0; -1 0 1];

% Mass matrix (\int w_j w_i dx)
matM = [ 1/12 1/24 1/24 ;
         1/24 1/12 1/24 ;
         1/24 1/24 1/12 ] * detJ;

% Stiffness matrix (\int Grad(w_j) Grad(w_i) dx)
matK = (gradW'*gradW) * detJ * 0.5;

% % Differentiation matrices (\int w_j d_x(w_i) ds)
matDX = gradW(1,:)' * ones(1,3) * (1/6)*detJ;
matDY = gradW(2,:)' * ones(1,3) * (1/6)*detJ;

end

% % Jacobian ([V1,V2,V3] <=> [(-1,-1),(1,-1),(-1,1)])
% J = [(V2-V1)' (V3-V1)']/2;
% detJ = abs(det(J));

% % Mass matrix (\int w_i w_j dx)
% matM = [ 1/3. 1/6. 1/6. ;
%         1/6. 1/3. 1/6. ;
%         1/6. 1/6. 1/3. ] * detJ;

% % Stiffness matrix (\int Grad(w_i) Grad(w_j) dx)
% V = [V2-V3; V3-V1; V1-V2];
% matK = 1/8 * (V*V') * 1/detJ;