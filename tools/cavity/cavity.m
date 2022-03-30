function [u,dx,dy] = cavity(x,y)
    global k;
    N=51;
    
    u   = zeros(size(x));
    dx  = zeros(size(x));
    dy  = zeros(size(x));
    d2x = zeros(size(x));
    d2y = zeros(size(x));
    
    for i=1:2:N
        for j=1:2:N
            ipi = i*pi;
            jpi = j*pi;
            
            ipi2 = ipi*ipi;
            jpi2 = jpi*jpi;
            
            c = 16./(ipi2+jpi2-k*k)/(ipi*jpi);
            
            u   = u    +      c*sin(ipi*x).*sin(jpi*y);
            
            dx  = dx   + ipi *c*cos(ipi*x).*sin(jpi*y);
            dy  = dy   + jpi *c*sin(ipi*x).*cos(jpi*y);
        end
    end
end
