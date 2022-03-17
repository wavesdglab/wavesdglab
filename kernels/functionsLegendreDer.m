% Assumption: x in [-1,1]

function val = functionsLegendreDer(x,degree)

legendre = functionsLegendre(x,degree);
x = x(:);

N = degree+1;
val = zeros(size(x(:),1), N);
val(:,1) = 0;  % order 0
for n=1:(N-1)
    val(:,n+1) = n*legendre(:,n) + x.*legendre(:,n);  % order n
end

end