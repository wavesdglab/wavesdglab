function mesh = setupBenchmark1D(xL, xR, numE)

if(~exist('benchmarks/1D','dir'))
    error('Error - No valid benchmark has been set.')
end

restoredefaultpath();
setup();
addpath('benchmarks/1D');
mesh = buildMesh1D(xL, xR, numE);

end