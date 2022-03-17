% Assumption: x,y in [-1,1]

function [valDX, valDY] = functionsShapeDerTRI(x,y,degree)

assert(size(x,1) == size(y,1));
L = size(x,1);

N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;

l1 =  (y+1)/2;
l2 = -(x+y)/2;
l3 =  (x+1)/2;

l1DX = 0   * ones(L,1);
l1DY = 0.5 * ones(L,1);
l2DX =-0.5 * ones(L,1);
l2DY =-0.5 * ones(L,1);
l3DX = 0.5 * ones(L,1);
l3DY = 0   * ones(L,1);

valDX = zeros(L,N);
valDY = zeros(L,N);
n = 1;

% node 1
valDX(:,n) = l2DX;
valDY(:,n) = l2DY;
n=n+1;

% node 2
valDX(:,n) = l3DX;
valDY(:,n) = l3DY;
n=n+1;

% node 3
valDX(:,n) = l1DX;
valDY(:,n) = l1DY;
n=n+1;

% edge 1
Nedg = degree+2; % to modify
xloc   = l3  -l2;
xlocDX = l3DX-l2DX;
xlocDY = l3DY-l2DY;
coef   = 4*l2  .*l3 ./ (1. - xloc.^2);
coefDX = 4*l2DX.*l3 ./ (1. - xloc.^2) + 4*l2.*l3DX ./ (1. - xloc.^2) - 4*l2.*l3 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDX);
coefDY = 4*l2DY.*l3 ./ (1. - xloc.^2) + 4*l2.*l3DY ./ (1. - xloc.^2) - 4*l2.*l3 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDY);
valEdg    = functionsLobbato(xloc,    degree);
valEdgDer = functionsLobbatoDer(xloc, degree);
size(valEdg)
for ne = 1:degree-1
    valDX(:,n) = coefDX .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDX;
    valDY(:,n) = coefDY .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDY;
    n=n+1;
end

% edge 2
Nedg = degree+2; % to modify
xloc   = l1  -l3;
xlocDX = l1DX-l3DX;
xlocDY = l1DY-l3DY;
coef   = 4*l3  .*l1 ./ (1. - xloc.^2);
coefDX = 4*l3DX.*l1 ./ (1. - xloc.^2) + 4*l3.*l1DX ./ (1. - xloc.^2) - 4*l3.*l1 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDX);
coefDY = 4*l3DY.*l1 ./ (1. - xloc.^2) + 4*l3.*l1DY ./ (1. - xloc.^2) - 4*l3.*l1 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDY);
valEdg    = functionsLobbato(xloc,    degree);
valEdgDer = functionsLobbatoDer(xloc, degree);
for ne = 1:degree-1
    valDX(:,n) = coefDX .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDX;
    valDY(:,n) = coefDY .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDY;
    n=n+1;
end

% edge 3
Nedg = degree+2; % to modify
xloc   = l2  -l1;
xlocDX = l2DX-l1DX;
xlocDY = l2DY-l1DY;
coef   = 4*l1  .*l2 ./ (1. - xloc.^2);
coefDX = 4*l1DX.*l2 ./ (1. - xloc.^2) + 4*l1.*l2DX ./ (1. - xloc.^2) - 4*l1.*l2 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDX);
coefDY = 4*l1DY.*l2 ./ (1. - xloc.^2) + 4*l1.*l2DY ./ (1. - xloc.^2) - 4*l1.*l2 ./ (1. - xloc.^2).^2 .* (- 2*xloc.*xlocDY);
valEdg    = functionsLobbato(xloc,    degree);
valEdgDer = functionsLobbatoDer(xloc, degree);
for ne = 1:degree-1
    valDX(:,n) = coefDX .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDX;
    valDY(:,n) = coefDY .* valEdg(:,2+ne) + coef .* valEdgDer(:,2+ne) .* xlocDY;
    n=n+1;
end

% face
for n1 = 1:degree-1
    for n2 = 1:degree-1-n1
        valDX(:,n) = l1DX .* l2.^n1 .* l3.^n2 + l1 .* (n1 .* l2.^(n1-1) .* l2DX) .* l3.^n2 + l1 .* l2.^n1 .* (n2 .* l3.^(n2-1) .* l3DX);
        valDY(:,n) = l1DY .* l2.^n1 .* l3.^n2 + l1 .* (n1 .* l2.^(n1-1) .* l2DY) .* l3.^n2 + l1 .* l2.^n1 .* (n2 .* l3.^(n2-1) .* l3DY);
        n=n+1;
    end
end

end