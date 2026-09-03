function setup()

addpath(genpath('src/'));
addpath(genpath('tools/'));
addpath('scripts/');

% directoryGmsh = ':/Applications/Gmsh.app/Contents/MacOS/';
% directoryGmsh = ':~/Desktop/gmsh-4.11.1-Linux64/bin';

path1 = getenv('PATH');
if(~contains(path1, directoryGmsh))
    path1 = [path1 directoryGmsh];
    setenv('PATH', path1);
end

setenv('PATH', path1);

set(0, 'DefaultLineLineWidth', 2);

end