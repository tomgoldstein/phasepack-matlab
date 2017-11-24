%                  runBenchmarkTransMatGDvsLBFGS.m
% 
% This example will compare the convergence of truncated Wirtinger flow
% with regular gradient descent to trunated Wirtinger flow with L-BFGS
% acceleration.  The measurement matrix used for the comparison is the
% empirical transmission matrix.  
%
% Note: you will need to download the transmission matrix datasets to run
% this script.  See the user guide for instructions.
% 
% This script does the following:
% 
% 1. Set up parameters and create a list of algorithm structs. 
%
% 2. Invoke the general benchmark function benchmarkPR. A graph of errors
% (under specified error metrics) of different algorithms at each number of
% iterations will be shown.
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017


%% -----------------------------START----------------------------------


clc
clear
close all
%% 1.Set up parameters
% Choose x label (values shown on the x axis of the benchmark plot) and 
% y label (values shown on the y-axis). The value on the x axis is the
% number of iterations.
xitem = 'iterations';
xvalues = [10 50 100 500 1000]; % Iterations allowed
yitem= 'reconerror';


% Choose Dataset
dataSet = 'transmissionMatrix';

% Set up general parameters
params.verbose = false;
params.numTrials = 20;        % run several random trials for each scenario, and report average results
params.n = 256;               % num of unknown elements
params.m = 20*params.n;       % number of measurements
params.isComplex = true;      % use complex matrices? or just stick to real?
params.policy = 'median';     % report median accuracy over the random trials


% Create two different versions of truncated wirtinger flow.  One using a
% steepest descent optimizer, and one using L-BFGS.  We also specify a
% 'label' for each algorithm, which is used to produce the legend of the
% plot.
twf_sd = struct('algorithm','twf','searchMethod','steepestDescent','label','TWF-SD');
twf_lbfgs = struct('algorithm','twf','searchMethod','LBFGS','label','TWF-LBFGS');                                       


% Grab your pick of algorithms.
algorithms = {twf_sd, twf_lbfgs};


% Run benchmark
benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params);
