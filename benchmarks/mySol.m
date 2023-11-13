function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k TAGbench R_disk L L_PML

switch TAGbench
    case 'cavity'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'waveguide'
        L = 1.;
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
    case 'scatteringPML'
        solU = solScattPlaneWaveHard(k,R_disk,x,y);
        solF = 0*x;
        solDx = -1i*k*exp(1i*k*x);
        solDy = 0*x;
        solVx = 0*x;
        solVy = 0*x;
    case 'scattering_rect'
        % no solU 
        theta = 4*pi/10;
        solU = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF = 0*x;
        solDx = -1i*k*cos(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solDy = -1i*k*sin(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solVx = 0*x;
        solVy = 0*x;
    case 'scattering_square'
        % no solU 
        theta = 4*pi/10;
        solU = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF = 0*x;
        solDx = -1i*k*cos(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solDy = -1i*k*sin(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solVx = 0*x;
        solVy = 0*x;
    otherwise
        warning('Error - No valid benchmark has been set.')
end
end