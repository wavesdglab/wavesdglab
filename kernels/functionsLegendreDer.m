% Assumption: x in [-1,1]

function val = functionsLegendreDer(x,N)

legendre = functionsLegendre(x,N);
x = x(:);

val = zeros(size(x(:),1), N);
val(:,1) = 0;  % order 0
for n=1:(N-1)
    val(:,n+1) = n*legendre(:,n) + x.*legendre(:,n);  % order n
end

end