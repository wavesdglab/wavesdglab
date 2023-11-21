% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function setup()

addpath(genpath('benchmarks/'));
addpath(genpath('scripts/'));
addpath(genpath('sources/'));
addpath(genpath('tools/'));
addpath(genpath('tools/eigtool/'));

path1 = getenv('PATH');
if(~contains(path1, ':/Applications/Gmsh.app/Contents/MacOS/'))
    path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
    setenv('PATH', path1);
end

set(0, 'DefaultLineLineWidth', 2);

end
