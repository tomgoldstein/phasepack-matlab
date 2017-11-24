%                        runBenchmarkTransMatMN.m
% 
% This example benchmarks algorithms based on their ability to reconstruct
% a synthetic signal (random Gaussian) using measurements acquired from a
% real data matrix (a transmission matrix).  The benchmark shows how the
% different methods behave as the number of rows sampled from the
% transmission matrix (i.e., the number 'm' of measurements) increases.
%
% This benchmark looks at the behavior of the methods using a *real*
% empirical measurement matrix, but using synthetic signals so that
% reconstruction accuracy can be measured directly.
%
% Note:  You will need to download the empirical datasets before this
% benchmark can be run.  See the user guide.
% 
% Expected time to run: 5mins
%
% The script does the following:
% 
% 1. Set up parameters and create a list of algorithm structs. 

% 2. Invoke the general benchmark function benchmarkPR. A graph of errors
% (under specified error metrics) of different algorithms will be shown as
% the number of measurements varies.
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
xitem = 'm/n';      % Vary the number of samples used for reconstruction
xvalues = [1 1.5 2 2.5 3 3.5 4 4.5 5 7.5 10 12.5 15]; % The over-sampling ratios to try. I.e, let m=3*n, then 5*n, etc...
yitem = 'reconerror'; % Plot the reconstruction error of each method

% Choose dataset: in this case we use the empirical transmission matrix
dataSet = 'transmissionMatrix';

% Set up general parameters
params.verbose = false;
params.numTrials = 10;        % run several random trials for each scenario, and report average results
params.n = 256;              % num of unknown elements.  When using the 
                             % transmissionMatrix dataset, this can be one of {246,1600,4096}
params.policy = 'median';    % report the median performance over the random trials


% Create a list of algorithms structs
% You can specify the algorithm name, the initializer, and other options
% for each method in the struct.
wf = struct('initMethod','spectral','algorithm','wirtflow');
twf = struct('algorithm','twf'); 
rwf = struct('algorithm','rwf');
ampflow = struct('algorithm','amplitudeflow');
taf = struct('initMethod','orthogonal','algorithm','taf');
raf = struct('initMethod','weighted','algorithm','raf');
fienup = struct('algorithm','fienup');
gs = struct('algorithm','gerchbergsaxton');
cd = struct('algorithm','coordinatedescent','maxIters',3000);
kac = struct('algorithm','kaczmarz','maxIters',1000);
pmax = struct( 'algorithm','phasemax', 'maxIters',1000);
plamp = struct('algorithm', 'phaselamp');
scgm = struct('algorithm','sketchycgm');     
plift = struct('algorithm','phaselift','maxIters',1000);                                             


% Grab your pick of algorithms to benchmark.
algorithms = {wf,raf,fienup,gs,pmax,plamp,plift};

% Run benchmark
benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params);
