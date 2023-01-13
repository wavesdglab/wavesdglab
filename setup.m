function setup()

addpath(genpath('benchmarks1D/'));
addpath(genpath('benchmarks2D/'));
addpath(genpath('kernels/'));
addpath(genpath('scripts/'));
addpath(genpath('tools/'));

path1 = getenv('PATH');
if(~contains(path1, ':/Applications/Gmsh.app/Contents/MacOS/'))
    path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
    setenv('PATH', path1);
end

set(0, 'DefaultLineLineWidth', 2);

end