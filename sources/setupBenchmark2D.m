function mesh = setupBenchmark2D(tag)

if(~exist(strcat('benchmarks/',tag),'dir'))
    error('Error - No valid benchmark has been set.')
end

restoredefaultpath();
setup();
addpath(strcat('benchmarks/',tag));
mesh = myBenchmark();

end