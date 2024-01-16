% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function setup()

addpath(genpath('benchmarks/'));
addpath(genpath('scripts/'));
addpath(genpath('sources/'));
addpath(genpath('tools/'));

directoryGmsh = ':/Applications/Gmsh.app/Contents/MacOS/';
% directoryGmsh = ':/home/pescuma/Desktop/gmsh-4.11.1-Linux64/bin';

if(exist('directoryGmsh'))
    path1 = getenv('PATH');
    if(~contains(path1, directoryGmsh))
        path1 = [path1 directoryGmsh];
        setenv('PATH', path1);
    end
else
    error('Variable "directoryGmsh" must be defined in "setup.m".')
end

set(0, 'DefaultLineLineWidth', 2);

end