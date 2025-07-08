% Assumption: x in [-1,1]

function val = functionsJacobiDer(x,alpha,beta,N)

x = x(:);
valJacobi = functionsJacobi(x,alpha,beta,N);

val = zeros(size(x,1),N);
val(:,1) = 0;
for n=1:(N-1)
    b1 = (2*n+alpha+beta)*(1-x.^2);
    b2 = n*(alpha-beta-(2*n+alpha+beta)*x);
    b3 = 2*(n+alpha)*(n+beta);
    val(:,n+1) = (b2.*valJacobi(:,n+1) + b3.*valJacobi(:,n)) ./ b1;
end

end