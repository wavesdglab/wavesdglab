function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M rho c omega;

souU = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

% % LEFT + RIGHT PROPAGATING WAVES
% souU = (1+M)/2*exp(1i*omega*x/(c*(1+M)))+(1-M)/2*exp(-1i*omega*x/(c*(1-M)));
% souDx = 1i*omega/(2*c) * (exp(1i*omega*x/(c*(1+M)))-exp(-1*omega*x/(c*(1-M))));
% souVx = 1/(2*rho*c)*((1+M)*exp(1i*omega*x/(c*(1+M)))-(1-M)*exp(-1i*omega*x/(c*(1-M))));

% ONLY RIGHT PROPAGATING WAVE
souU = (1+M)/2*exp(1i*omega*x/(c*(1+M)));
souDx = 1i*omega/(2*c) * exp(1i*omega*x/(c*(1+M)));
souVx = (1+M)/(2*rho*c) * exp(1i*omega*x/(c*(1+M)));

end