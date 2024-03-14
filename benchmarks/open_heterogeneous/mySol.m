function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global eta1 eta2 k1 k2 c1 c2

thetaI = pi/6;                           % thetaI \in [0, \pi/2]
sinI = sin(thetaI);
cosI = cos(thetaI);
sinT = c2 / c1 * sinI;
cosT = sqrt(1-sinT^2);

if c2/c1 * sinI > 1
    error('Total reflection!')
end

I = 1;
R = I * (eta2*cosI - eta1*cosT) / (eta1 * cosT + eta2 * cosI) * exp(1i*k1*cosI);
T = I * (2*eta2*cosI) / (eta1*cosT + eta2*cosI) * exp(1i*(k1*cosI - k2*cosT)/2);

if (min(x)<0.5)
    solU = I * exp(1i*k1*(cosI*x+sinI*y)) + R * exp(1i*k1*(-cosI*x+sinI*y));
    solDx = 1i * k1 * cosI * (I * exp(1i*k1*(cosI*x+sinI*y)) - R * exp(1i*k1*(-cosI*x+sinI*y)));
    solDy = 1i * k1 * sinI * (I * exp(1i*k1*(cosI*x+sinI*y)) + R * exp(1i*k1*(-cosI*x+sinI*y)));
    solF = 0*x + 0*y;
    solVx = 1/(1i*k1*eta1) * solDx;
    solVy = 1/(1i*k1*eta1) * solDy;
else
    solU = T * exp(1i*k2*(cosT*x + sinT*y));
    solDx = 1i * k2 * cosT * T * exp(1i*k2*(cosT*x + sinT*y));
    solDy = 1i * k2 * sinT * T * exp(1i*k2*(cosT*x + sinT*y));
    solF = 0*x + 0*y;
    solVx = 1/(1i*k2*eta2) * solDx;
    solVy = 1/(1i*k2*eta2) * solDy;
end

end