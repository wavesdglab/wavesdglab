% Assumption: x in [-1,1]

function val = functionsLobbatoDer(x,N)

% val = functionsBernsteinDer(x,N);

legendre = functionsLegendre(x,N);
x = x(:);

val = zeros(size(x(:),1),N);
val(:,1) = -1/2;  % order 0
val(:,2) = +1/2;  % order 1
for n=2:(N-1)
    val(:,n+1) = sqrt(n - 1/2) * legendre(:,n);  % order n
end

end