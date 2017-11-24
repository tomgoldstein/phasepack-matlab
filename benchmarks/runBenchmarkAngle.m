%                           runBenchmarkAngle.m
% 
% This example benchmarks algorithms based on their ability to reconstruct
% a synthetic signal (random Gaussian) using synthetic measurements 
% (random Gaussian).  The benchmark shows how the
% different methods behave as the angle between the initalizer and the
% true signal changes.
% 
% This script does the following:
% 
% 1. Set up parameters and create a list of algorithm structs. 
%
% 2. Invoke the general benchmark function benchmarkPR. A graph of errors
% (under specified error metrics) of different algorithms will be shown.
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
xitem = 'angle';
xvalues = 0.001:pi/10:pi/2; 
yitem = 'reconerror';


% Choose Dataset and set up dataSet '1DGaussian' specific parameters
dataSet = '1DGaussian';

% Set up general parameters
params.verbose = false;
params.numTrials = 10;       % run several random trials for each scenario, and report average results
params.n = 100;              % num of unknown elements
params.m = 5*params.n;      % num of measurements
params.isComplex = true;    % use complex matrices? or just stick to real?
params.policy = 'median';	% use the median performance over trails to assess algortihms


% Create a list of algorithms structs
wf = struct('algorithm','wirtflow');
twf = struct('algorithm','twf'); 
rwf = struct('algorithm','rwf');
ampflow = struct('algorithm','amplitudeflow');
taf = struct('algorithm','taf');
raf = struct('algorithm','raf');
fienup = struct('algorithm','fienup');
gs = struct('algorithm','gerchbergsaxton');
cd = struct('algorithm','coordinatedescent','maxIters',3000);
kac = struct('algorithm','kaczmarz','maxIters',1000);
pmax = struct( 'algorithm','phasemax', 'maxIters',1000);
plamp = struct('algorithm', 'phaselamp');
scgm = struct('algorithm','sketchycgm');     
plift = struct('algorithm','phaselift','maxIters',1000);                                             

% Grab your pick of algorithms.
algorithms = {rwf, pmax,plamp};

% Run benchmark
benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params);
