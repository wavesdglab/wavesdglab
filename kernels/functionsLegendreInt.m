% Assumption: x in [-1,1]

function val = functionsLegendreInt(x,N)

legendre = functionsLegendre(x,N);
x = x(:);

val = zeros(size(x(:),1),N);
val(:,1) = 0;  % order 0
for n=1:(N-1)
    val(:,n+1) = (x.*legendre(:,n+1) - legendre(:,n)) / (n+1);  % order n
end

end