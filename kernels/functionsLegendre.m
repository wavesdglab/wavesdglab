% Assumption: x in [-1,1]

function val = functionsLegendre(x,N)

x = x(:);

val = zeros(size(x(:),1),N);
val(:,1) = 1;  % order 0
val(:,2) = x;  % order 1
for n=2:(N-1)
    val(:,n+1) = ((2*n-1)*x.*val(:,n) - (n-1)*val(:,n-1)) / n;  % order n
end

end