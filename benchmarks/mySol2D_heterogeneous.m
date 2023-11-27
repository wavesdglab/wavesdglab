function [solU, solDx, solDy, solF, solVx, solVy, rho, c, eta] = mySol2D_heterogeneous(x,y)

global TAGbench omega;

% Physical parameters
rho1 = 1;
c1 = 2;
rho2 = 1;
c2 = 0.8;

eta1 = rho1 * c1;
eta2 = rho2 * c2;
k1 = omega / c1;
k2 = omega / c2;

if min(x) < 0.5
    rho = rho1;
    c = c1;
else
    rho = rho2;
    c = c2;
end
eta = rho * c;
% k = omega / c;

switch TAGbench
    case 'cavity'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'waveguide'
        L = 4.;
        theta = 30.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'open'
        theta = pi/4;
        solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF  = 0*x;
        solDx = 1i*k*cos(theta) * solU;
        solDy = 1i*k*sin(theta) * solU;
        solVx = cos(theta) * solU;
        solVy = sin(theta) * solU;
    case 'cavity (heterogeneous)'

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

    case 'waveguide (heterogeneous)'
        L = 4.;
        theta = 45.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*eta*k);
        solVy = solDy/(1i*eta*k);
    case 'open (heterogeneous)'

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

    otherwise
        warning('Error - No valid benchmark has been set.')
end

end