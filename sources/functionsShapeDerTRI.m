% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x,y in [-1,1]

function [valDX, valDY] = functionsShapeDerTRI(x,y,degree)

assert(size(x,1) == size(y,1));

L = size(x,1);
N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;

% Barycentric coordinates
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

% Nodal modes
valDX(:,n) = l2DX; valDY(:,n) = l2DY; n=n+1;
valDX(:,n) = l3DX; valDY(:,n) = l3DY; n=n+1;
valDX(:,n) = l1DX; valDY(:,n) = l1DY; n=n+1;

% Edge modes
kernel1 = functionsJacobi(l3-l2,1,1,degree-1);
kernel2 = functionsJacobi(l1-l3,1,1,degree-1);
kernel3 = functionsJacobi(l2-l1,1,1,degree-1);
kernelDer1 = functionsJacobiDer(l3-l2,1,1,degree-1);
kernelDer2 = functionsJacobiDer(l1-l3,1,1,degree-1);
kernelDer3 = functionsJacobiDer(l2-l1,1,1,degree-1);
for ne = 1:degree-1
    valDX(:,n) = (l2DX .* l3 + l2 .* l3DX) .* kernel1(:,ne) ...
        + l2 .* l3 .* kernelDer1(:,ne) .* (l3DX-l2DX);
    valDY(:,n) = (l2DY .* l3 + l2 .* l3DY) .* kernel1(:,ne) ...
        + l2 .* l3 .* kernelDer1(:,ne) .* (l3DY-l2DY);
    n=n+1;
end
for ne = 1:degree-1
    valDX(:,n) = (l3DX .* l1 + l3 .* l1DX) .* kernel2(:,ne) ...
        + l3 .* l1 .* kernelDer2(:,ne) .* (l1DX-l3DX);
    valDY(:,n) = (l3DY .* l1 + l3 .* l1DY) .* kernel2(:,ne) ...
        + l3 .* l1 .* kernelDer2(:,ne) .* (l1DY-l3DY);
    n=n+1;
end
for ne = 1:degree-1
    valDX(:,n) = (l1DX .* l2 + l1 .* l2DX) .* kernel3(:,ne) ...
        + l1 .* l2 .* kernelDer3(:,ne) .* (l2DX-l1DX);
    valDY(:,n) = (l1DY .* l2 + l1 .* l2DY) .* kernel3(:,ne) ...
        + l1 .* l2 .* kernelDer3(:,ne) .* (l2DY-l1DY);
    n=n+1;
end

% Face modes
for n1 = 1:degree-1
    for n2 = 1:degree-1-n1
        valDX(:,n) = (l1DX.*l2.*l3 + l1.*l2DX.*l3 + l1.*l2.*l3DX) .* kernel1(:,n1) .* kernel3(:,n2) ...
            + l1.*l2.*l3 .* kernelDer1(:,n1) .* kernel3(:,n2) .* (l3DX-l2DX) ...
            + l1.*l2.*l3 .* kernel1(:,n1) .* kernelDer3(:,n2) .* (l2DX-l1DX);
        valDY(:,n) = (l1DY.*l2.*l3 + l1.*l2DY.*l3 + l1.*l2.*l3DY) .* kernel1(:,n1) .* kernel3(:,n2) ...
            + l1.*l2.*l3 .* kernelDer1(:,n1) .* kernel3(:,n2) .* (l3DY-l2DY) ...
            + l1.*l2.*l3 .* kernel1(:,n1) .* kernelDer3(:,n2) .* (l2DY-l1DY);
        n=n+1;
    end
end

end