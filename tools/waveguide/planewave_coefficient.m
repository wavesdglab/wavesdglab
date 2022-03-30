function g = planewave_coefficient(n)
	threshold = 1.e-10;

	global k;
	global theta;

	npi = n*pi;

	a = k*sin(theta);

	ip = i*(a+npi);
	im = i*(a-npi);

	Ip = 1.;
	Im = 1.;

	if (abs(ip) > threshold)
		Ip = (exp(ip)-1.)/ip;
	end

	if (abs(im) > threshold)
		Im = (exp(im)-1.)/im;
	end

	g = (Ip-Im)/(2.*i);
end
