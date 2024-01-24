% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function setup()

addpath(genpath('benchmarks/'));
addpath(genpath('scripts/'));
addpath(genpath('sources/'));
addpath(genpath('tools/'));

path1 = getenv('PATH');
if(~contains(path1, ':/Applications/Gmsh.app/Contents/MacOS/'))
    path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
    setenv('PATH', path1);
end

% path1 = getenv('PATH');
% path1 = [':/home/pescuma/Desktop/gmsh-4.11.1-Linux64/bin' path1];

setenv('PATH', path1);

set(0, 'DefaultLineLineWidth', 2);

end