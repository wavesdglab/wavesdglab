% Assumption: x in [-1,1]

function val = functionsLobbatoDer(x,degree)

legendre = functionsLegendre(x,degree);
x = x(:);

N = degree+1;
val = zeros(size(x(:),1),N);

% nodal modes
val(:,1) = -1/2;
val(:,2) = +1/2;

% edge modes
for n=2:(N-1)
    val(:,n+1) = sqrt(n - 1/2) * legendre(:,n);
end

end