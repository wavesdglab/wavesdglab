% Assumption: x in [-1,1]

function val = functionsLagrangeTRI(x,y,degree)

assert(size(x,1) == size(y,1));

L = size(x,1);
N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;

% Barycentric coordinates
l1 = -(x+y)/2;
l2 =  (x+1)/2;
l3 =  (y+1)/2;

val = zeros(L,N);
n = 1;

% P1
if(degree == 1)
    val(:,n) = l1; n=n+1;
    val(:,n) = l2; n=n+1;
    val(:,n) = l3; n=n+1;
end

% P2
if(degree == 2)
    val(:,n) = l1.*(2*l1-1); n=n+1;
    val(:,n) = l2.*(2*l2-1); n=n+1;
    val(:,n) = l3.*(2*l3-1); n=n+1;
    val(:,n) = 4*l1.*l2; n=n+1;
    val(:,n) = 4*l2.*l3; n=n+1;
    val(:,n) = 4*l1.*l3; n=n+1;
end

% P3
if(degree == 3)
    val(:,n) = 0.5*l1.*(3*l1-1).*(3*l1-2); n=n+1;
    val(:,n) = 0.5*l2.*(3*l2-1).*(3*l2-2); n=n+1;
    val(:,n) = 0.5*l3.*(3*l3-1).*(3*l3-2); n=n+1;
    val(:,n) = (9/2)*l1.*(3*l1-1).*l2; n=n+1;
    val(:,n) = (9/2)*l2.*(3*l2-1).*l1; n=n+1;
    val(:,n) = (9/2)*l2.*(3*l2-1).*l3; n=n+1;
    val(:,n) = (9/2)*l3.*(3*l3-1).*l2; n=n+1;
    val(:,n) = (9/2)*l3.*(3*l3-1).*l1; n=n+1;
    val(:,n) = (9/2)*l1.*(3*l1-1).*l3; n=n+1;
    val(:,n) = 27*l1.*l2.*l3; n=n+1;
end

if(degree > 3)
    error('Error: functionsLagrangeTRI() not available for degree P=%i', degree);
end

end