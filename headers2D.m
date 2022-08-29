function headers2D()

addpath('benchmarks2D/');
addpath('kernels/');
addpath('tools/');
addpath('tools/triasymq/');
addpath('tools/cavity/');
addpath('tools/waveguide/');

path1 = getenv('PATH');
path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
setenv('PATH', path1);

set(0, 'DefaultLineLineWidth', 2);

end
