function [solU, solDx, solDy, solF, solVx, solVy] = mySol_heterogeneous(x,y,k)

global TAGbench omega;

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
    case 'cavity_heterogeneous'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'waveguide_heterogeneous'
        L = 4.;
        theta = 45.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'open_heterogeneous'
        
        rho1 = 1;
        c1 = 1;
        rho2 = 1;
        c2 = 1;

        thetaI = pi/4;                           % thetaI \in [0, \pi/2]
        sinI = sin(thetaI);
        cosI = cos(thetaI);
        sinT = c2 / c1 * sinI;
        cosT = sqrt(1-sinT^2);

        
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;

        I = 1;
        R = I * (eta2*cosI - eta1*cosT) / (eta1 * cosT + eta2 * cosI) * exp(1i*k1*cosI);
        T = I * (2*eta2*cosI) / (eta1*cosT + eta2*cosI) * exp(1i*(k1*cosI - k2*cosT)/2);

        if (min(x)<1/2)
            solU = I * exp(1i*k1*(cosI*x + sinI*y)) + R * exp(1i*k1*(-cosI*x+sinI*y));
            solDx = 1i * k1 * cosI * (I * exp(1i*k1*(cosI*x + sinI*y)) - R * exp(1i*k1*(-cosI*x+sinI*y)));
            solDy = 1i * k1 * sinI * (I * exp(1i*k1*(cosI*x + sinI*y)) + R * exp(1i*k1*(-cosI*x+sinI*y)));
            solF = 0*x + 0*y;
            solVx = 1/(1*k1*eta1) * solDx;
            solVy = 1/(1*k1*eta1) * solDy;
        else
            solU = T * exp(1i*k2*(cosT*x + sinT*y));
            solDx = 1i * k2 * cosT * T * exp(1i*k2*(cosT*x + sinT*y));
            solDy = 1i * k2 * sinT * T * exp(1i*k2*(cosT*x + sinT*y));
            solF = 0*x + 0*y;
            solVx = 1/(1*k1*eta1) * solDx;
            solVy = 1/(1*k1*eta1) * solDy;
        end




%         theta = pi/4;
%         solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
%         solF  = 0*x;
%         solDx = 1i*k*cos(theta) * solU;
%         solDy = 1i*k*sin(theta) * solU;
%         solVx = cos(theta) * solU;
%         solVy = sin(theta) * solU;
    otherwise
        warning('Error - No valid benchmark has been set.')
end

end