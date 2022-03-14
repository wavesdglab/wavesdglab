% Build elemental matrices on the LINE

function [matM, matK, matD] = buildMatrixElemLIN(V1, V2, degree)

N = degree+1;
Q = ceil((2*(N-1)+1)/2);
[nodes, weights] = quadratureGaussLIN(Q);
valQ = functionsShapeLIN(nodes,N);
derQ = functionsShapeDerLIN(nodes,N);

% Reference line [-1,1]
matMref = zeros(N,N);
matKref = zeros(N,N);
matDref = zeros(N,N);
for i=1:N
    for j=1:N
        matMref(i,j) = weights' * (valQ(:,i) .* valQ(:,j));
        matKref(i,j) = weights' * (derQ(:,i) .* derQ(:,j));
        matDref(i,j) = weights' * (valQ(:,i) .* derQ(:,j));
    end
end

matMref = [ 2/6 1/6 ; 1/6 2/6 ];
matKref = [ 1 -1 ; -1 1 ];
matDref = [ -1/2 -1/2 ; 1/2 1/2 ];

% Physical line [V1,V2]
J = norm(V2-V1);
matM = matMref*J;
matK = matKref/J;
matD = matDref;

end

% Definition of the elements of the matrices:
% - Mass matrix: \int w_i w_j ds
% - Stiffness matrix: \int d_s(w_i) d_s(w_j) ds
% - Differentiation matrix: \int w_i d_s(w_j) ds

% Particular case (N=2) :
% - matM = [ 2/6 1/6 ; 1/6 2/6 ];
% - matK = [ 1 -1 ; -1 1 ];
% - matD = [ -1/2 -1/2 ; 1/2 1/2 ];

% Mapping with the physical LINE [V1,V2] using Jacobian J = norm(V2-V1)
% - matMloc = matM*J;
% - matKloc = matK/J;
% - matDloc = matD;