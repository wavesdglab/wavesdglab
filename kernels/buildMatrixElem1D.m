function [matM, matK, matD] = buildMatrixElem1D(degree)

N = degree+1;
Q = ceil((2*(N-1)+1)/2);
[nodes, weights] = quadratureGauss1D(Q);

valQ = functionsShape1D(nodes,degree);
derQ = functionsShapeDer1D(nodes,degree);

matM = zeros(N,N);
matK = zeros(N,N);
matD = zeros(N,N);

for i=1:N
    for j=1:N
        matM(i,j) = weights' * (valQ(:,i) .* valQ(:,j));
        matK(i,j) = weights' * (derQ(:,i) .* derQ(:,j));
        matD(i,j) = weights' * (valQ(:,i) .* derQ(:,j));
    end
end

end