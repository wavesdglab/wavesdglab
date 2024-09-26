function [dsR_dt] = mySol_Robin(x,y,nx,ny)

global eta1 eta2 k1 k2 c1 c2

thetaI = pi/4;                           % thetaI \in [0, \pi/2]
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

    if (abs(nx) < 1e-5 && ny < 0 && abs(ny+1) < 1e-5)   % nx == 0 && ny == -1
        A = solDx;
        B = 0;
        C = 1i * k1 * sinI * solDx;
    end
    if (abs(nx) < 1e-5 && ny > 0 && abs(ny-1) < 1e-5)   % nx == 0 && ny == 1
        A = -solDx;
        B = 0;
        C = -1i * k1 * sinI * solDx;
    end
    if (nx < 0 && abs(nx+1) < 1e-5 && abs(ny) < 1e-5)   % nx == -1 && ny == 0
        A = -solDy;
        B = -1i * k1 * sinI * solDx;
        C = 0;
    end

    dsR_dt = A - 1/(1i * k1) * (nx * B + ny * C);

else
    solU = T * exp(1i*k2*(cosT*x + sinT*y));
    solDx = 1i * k2 * cosT * T * exp(1i*k2*(cosT*x + sinT*y));
    solDy = 1i * k2 * sinT * T * exp(1i*k2*(cosT*x + sinT*y));

    if (abs(nx) < 1e-5 && ny < 0 && abs(ny+1) < 1e-5)   % nx == 0 && ny == -1
        A = solDx;
        B = 0;
        C = 1i * k2 * sinI * solDx;
    end
    if (abs(nx) < 1e-5 && ny > 0 && abs(ny-1) < 1e-5)   % nx == 0 && ny == 1
        A = -solDx;
        B = 0;
        C = -1i * k2 * sinI * solDx;
    end
    if (nx > 0 && abs(nx-1) < 1e-5 && abs(ny) < 1e-5)   % nx == 1 && ny == 0
        A = solDy;
        B = 1i * k2 * sinI * solDx;
        C = 0;
    end

    dsR_dt = A - 1/(1i * k2) * (nx * B + ny * C);

end

end