%                           runBenchmarkSuccessRate.m
% 
% This example benchmarks algorithms based on their ability to reconstruct
% a synthetic signal (random Gaussian) using synthetic measurements 
% (random Gaussian).  The benchmark shows how the
% different methods behave as the number of measurements increases.  The
% y-axis plots the success rate (rate of exact signal recovery).  This is
% done by setting params.policy='successrate'.
%
% The algorithms currently used in this implementation are all instances of
% PhaseMax, but with different levels of accuracy in the initializer.  To
% control the level of initialization accuracy, we set the initializer to
% "angle", and specify the angle between the true signal, and the initial
% guess.
% 
% This script does the following:
% 
% 1. Set up parameters and create a list of algorithm structs. 
%
% 2. Invoke the general benchmark function benchmarkPR. A graph of errors
% (under specified error metrics) of different algorithms will be shown.
%
% The benchmark program compares the performance of specified algorithms on
% 1D gaussian measurements that have different m/n ratio(i.e. the ratio
% between the number of measurements and the number of unknowns).
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
% y label (values shown on the y-axis). The value on the x axis is the m/n
% ratio. The value on the y axis is 'reconerror', which is the relative
% 2-norm difference between the true and recovered signal.
xitem = 'm/n';
xvalues = 2:.5:12; 
yitem = 'reconerror';

% Choose Dataset and set up dataSet '1DGaussian' specific parameters
dataSet = '1DGaussian';

% Set up general parameters
params.verbose = false;
params.numTrials = 20;         % run several random trials for each scenario, and report average results
params.n = 500;                 % num of unknown elements
params.isComplex = true;        % use complex matrices? or just stick to real?
params.policy = 'successrate';	% use the successrate
params.successConstant = 1e-4;

%  Each of these algorithm is an instance of PhaseMax.  However, they each 
%  are initialized with starting points of different accuracies.  The 
%  "angle" initializer grabs the "initAngle" entry from the options, and
%  produces an initializer that makes this angle with the true signal.
pmax25 = struct('algorithm','phasemax','initMethod','angle');                                             
pmax25.tol=1e-6;
pmax25.initAngle=25/360*2*pi;

pmax36 = struct('algorithm','phasemax','initMethod','angle');                                             
pmax36.tol=1e-6;
pmax36.initAngle=36/360*2*pi;

pmax45 = struct('algorithm','phasemax','initMethod','angle');                                             
pmax45.tol=1e-6;
pmax45.initAngle=45/360*2*pi;

% Grab your pick of algorithms.
algorithms = {pmax25,pmax36,pmax45};


% Run benchmark
results = benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params);
