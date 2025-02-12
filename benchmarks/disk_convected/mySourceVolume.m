function val = mySourceVolume(x,y)

global M

xs=-M;
ys=0;

val = 0*x + 0*y;

if (xs>min(x) && xs<max(x) && ys>min(y) && ys<max(y))
    val = 0*x + 0*y + 1;
end