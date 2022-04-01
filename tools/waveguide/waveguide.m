function [u,dx,dy] = waveguide(x,y)
	N = 48;

	global k;
	global L;
	global theta;
    L = 4;
    theta = pi/3;

	u  = zeros(size(x));
	dx = zeros(size(x));
	dy = zeros(size(x));

	for n=1:N
		npi = n*pi;
		g = planewave_coefficient(n);

		[v,d] = ode(n,g,x);

		u  = u  +     v.*sin(npi*y);
		dx = dx +     d.*sin(npi*y);
		dy = dy + npi*v.*cos(npi*y);
	end
end
