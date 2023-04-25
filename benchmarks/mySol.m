function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k TAGbench R;

switch TAGbench
    case 'cavity'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'cavity2'
        [solU, solDx, solDy] = cavity2(x,y,k);
        % solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
        solF = -(x - 0.35) .* (y - 0.25) .* (40000 * ((x - 0.35).^2 + (y - 0.25).^2) + 1200 + k) .* exp(-((x - 0.35).^2 + (y - 0.25).^2) / (0.1^2));
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
    case 'scatteringHard'
        solU = solScattPlaneWaveHard(k,R,x,y);
        solF = exp(1i*k*x);
        % solF = 0*x;
        solDx = -1i*k*exp(1i*k*x);
        solDy = 0*x;
        % solVx = 0*x;
        % solVy = 0*x;
    case 'scatteringSoft'
        solU = solScattPlaneWaveSoft(k,R,x,y);
        solF = exp(1i*k*x);
        solDx = 0*x;
        solDy = 0*x;
        solVx = 0*x;
        solVy = 0*x;
    case 'scatteringPML_01'
        solU = solScattPlaneWaveHard(k,R,x,y);
        solF = exp(1i*k*x);
        solDx = -1i*k*exp(1i*k*x);
        solDy = 0*x;
        solVx = 0*x;
        solVy = 0*x;
    case 'scatteringPML_02'
        solU = solScattPlaneWaveHard(k,R,x,y);
        solF = exp(1i*k*x);
        solDx = -1i*k*exp(1i*k*x);
        solDy = 0*x;
    case 'scatteringPML_05'
        solU = solScattPlaneWaveHard(k,R,x,y);
        % solF = exp(1i*k*x);
        solF = 0*x;
        solDx = -1i*k*exp(1i*k*x);
        solDy = 0*x;
    otherwise
        warning('Error - No valid benchmark has been set.')
end
end