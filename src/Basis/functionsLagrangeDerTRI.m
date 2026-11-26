% Assumption: x in [-1,1]

function [valDX, valDY] = functionsLagrangeDerTRI(x,y,degree)

assert(size(x,1) == size(y,1));

L = size(x,1);
N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;

% Barycentric coordinates
l1 = -(x+y)/2;
l2 =  (x+1)/2;
l3 =  (y+1)/2;
l1DX =-0.5 * ones(L,1);
l1DY =-0.5 * ones(L,1);
l2DX = 0.5 * ones(L,1);
l2DY = 0   * ones(L,1);
l3DX = 0   * ones(L,1);
l3DY = 0.5 * ones(L,1);

valDX = zeros(L,N);
valDY = zeros(L,N);

% P1
if(degree == 1)
    valDX(:,1) = l1DX;
    valDX(:,2) = l2DX;
    valDX(:,3) = l3DX;
    valDY(:,1) = l1DY;
    valDY(:,2) = l2DY;
    valDY(:,3) = l3DY;
end

% P2
if(degree == 2)
    valDX(:,1) = l1DX.*(2*l1-1) + l1.*(2*l1DX);
    valDX(:,2) = l2DX.*(2*l2-1) + l2.*(2*l2DX);
    valDX(:,3) = l3DX.*(2*l3-1) + l3.*(2*l3DX);
    valDX(:,4) = 4*l1DX.*l2 + 4*l1.*l2DX;
    valDX(:,5) = 4*l2DX.*l3 + 4*l2.*l3DX;
    valDX(:,6) = 4*l1DX.*l3 + 4*l1.*l3DX;
    valDY(:,1) = l1DY.*(2*l1-1) + l1.*(2*l1DY);
    valDY(:,2) = l2DY.*(2*l2-1) + l2.*(2*l2DY);
    valDY(:,3) = l3DY.*(2*l3-1) + l3.*(2*l3DY);
    valDY(:,4) = 4*l1DY.*l2 + 4*l1.*l2DY;
    valDY(:,5) = 4*l2DY.*l3 + 4*l2.*l3DY;
    valDY(:,6) = 4*l1DY.*l3 + 4*l1.*l3DY;
end

% P3
if(degree == 3)
    valDX(:,1) = 0.5*l1DX.*(3*l1-1).*(3*l1-2) + 0.5*l1.*(3*l1DX).*(3*l1-2) + 0.5*l1.*(3*l1-1).*(3*l1DX);
    valDX(:,2) = 0.5*l2DX.*(3*l2-1).*(3*l2-2) + 0.5*l2.*(3*l2DX).*(3*l2-2) + 0.5*l2.*(3*l2-1).*(3*l2DX);
    valDX(:,3) = 0.5*l3DX.*(3*l3-1).*(3*l3-2) + 0.5*l3.*(3*l3DX).*(3*l3-2) + 0.5*l3.*(3*l3-1).*(3*l3DX);
    valDX(:,4) = (9/2)*l1DX.*(3*l1-1).*l2 + (9/2)*l1.*(3*l1DX).*l2 + (9/2)*l1.*(3*l1-1).*l2DX;
    valDX(:,5) = (9/2)*l1DX.*(3*l1-1).*l3 + (9/2)*l1.*(3*l1DX).*l3 + (9/2)*l1.*(3*l1-1).*l3DX;
    valDX(:,6) = (9/2)*l2DX.*(3*l2-1).*l1 + (9/2)*l2.*(3*l2DX).*l1 + (9/2)*l2.*(3*l2-1).*l1DX;
    valDX(:,7) = (9/2)*l2DX.*(3*l2-1).*l3 + (9/2)*l2.*(3*l2DX).*l3 + (9/2)*l2.*(3*l2-1).*l3DX;
    valDX(:,8) = (9/2)*l3DX.*(3*l3-1).*l1 + (9/2)*l3.*(3*l3DX).*l1 + (9/2)*l3.*(3*l3-1).*l1DX;
    valDX(:,9) = (9/2)*l3DX.*(3*l3-1).*l2 + (9/2)*l3.*(3*l3DX).*l2 + (9/2)*l3.*(3*l3-1).*l2DX;
    valDX(:,10) = 27*l1DX.*l2.*l3 + 27*l1.*l2DX.*l3 + 27*l1.*l2.*l3DX;
    valDY(:,1) = 0.5*l1DY.*(3*l1-1).*(3*l1-2) + 0.5*l1.*(3*l1DY).*(3*l1-2) + 0.5*l1.*(3*l1-1).*(3*l1DY);
    valDY(:,2) = 0.5*l2DY.*(3*l2-1).*(3*l2-2) + 0.5*l2.*(3*l2DY).*(3*l2-2) + 0.5*l2.*(3*l2-1).*(3*l2DY);
    valDY(:,3) = 0.5*l3DY.*(3*l3-1).*(3*l3-2) + 0.5*l3.*(3*l3DY).*(3*l3-2) + 0.5*l3.*(3*l3-1).*(3*l3DY);
    valDY(:,4) = (9/2)*l1DY.*(3*l1-1).*l2 + (9/2)*l1.*(3*l1DY).*l2 + (9/2)*l1.*(3*l1-1).*l2DY;
    valDY(:,5) = (9/2)*l1DY.*(3*l1-1).*l3 + (9/2)*l1.*(3*l1DY).*l3 + (9/2)*l1.*(3*l1-1).*l3DY;
    valDY(:,6) = (9/2)*l2DY.*(3*l2-1).*l1 + (9/2)*l2.*(3*l2DY).*l1 + (9/2)*l2.*(3*l2-1).*l1DY;
    valDY(:,7) = (9/2)*l2DY.*(3*l2-1).*l3 + (9/2)*l2.*(3*l2DY).*l3 + (9/2)*l2.*(3*l2-1).*l3DY;
    valDY(:,8) = (9/2)*l3DY.*(3*l3-1).*l1 + (9/2)*l3.*(3*l3DY).*l1 + (9/2)*l3.*(3*l3-1).*l1DY;
    valDY(:,9) = (9/2)*l3DY.*(3*l3-1).*l2 + (9/2)*l3.*(3*l3DY).*l2 + (9/2)*l3.*(3*l3-1).*l2DY;
    valDY(:,10) = 27*l1DY.*l2.*l3 + 27*l1.*l2DY.*l3 + 27*l1.*l2.*l3DY;
end

end