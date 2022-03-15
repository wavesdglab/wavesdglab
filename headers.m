function headers()

addpath('kernels/');
addpath('tools/');

path1 = getenv('PATH');
path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
setenv('PATH', path1);

end