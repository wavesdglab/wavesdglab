% Assumption: x in [-1,1]

function val = functionsLegendreInt(x,degree)

legendre = functionsLegendre(x,degree);
x = x(:);

N = degree+1;
val = zeros(size(x(:),1),N);
val(:,1) = 0;  % order 0
for n=1:(N-1)
    val(:,n+1) = (x.*legendre(:,n+1) - legendre(:,n)) / (n+1);  % order n
end

end