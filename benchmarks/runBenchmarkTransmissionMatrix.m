%                    runBenchmarkTransmissionMatrix.m
% 
% Benchmark algorthms using real/empirical data.  The measurement matrix is
% a transmission matrix, and the measurements were acquired using an
% optical aparatus. The data acquisition is described in "Coherent
% Inverse Scattering via Transmission Matrices: Efficient Phase Retrieval
% Algorithms and a Public Dataset." 
%    This script will reconstruct images using various algorithms, and
% report the image quality of the results.
%
% See the header of transmissionMatrixExperiment.m, or the paper cited
% above for more details.
%
% This script does the following: 
% 1. Set up parameters and create a list of
% algorithm structs. 
% 2. Invoke the general benchmark function benchmarkTransmissionMatrix.
%
% The benchmark program compares the performance of specified algorithms on
% reconstructing an image.  By default, this script reconstructs a 16x16
% image.  However, it can also reconstruct a 40x40 or 64x64 image by
% changing the options below.  Note, the runtime and memory requirements 
% are very high when reconstructing larger images.
% 
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


%clc
%clear
%close all

%% 1.Set up parameters
imageSize = 64;         % Side-length of image to reconstruct. Valid choices: 16,40,64
datasetSelection = 5;   % Which dataset to select measurements from.  Valid choices: 1-5
residualConstant = 0.3; 


% Create a list of algorithms structs
wf = struct('initMethod','spectral','algorithm','wirtflow');
twf = struct('algorithm','twf'); 
rwf = struct('algorithm','rwf');
ampflow = struct('algorithm','amplitudeflow');
taf = struct('initMethod','orthogonal','algorithm','taf');
raf = struct('initMethod','weighted','algorithm','raf');
fienup = struct('algorithm','fienup');
gs = struct('algorithm','gerchbergsaxton');
cd = struct('algorithm','coordinatedescent','maxIters',50000);
kac = struct('algorithm','kaczmarz','maxIters',1000);
pmax = struct( 'algorithm','phasemax','maxIters',1000);
plamp = struct('algorithm', 'phaselamp','maxIters',1000);
scgm = struct('algorithm','sketchycgm');     
plift = struct('algorithm','phaselift');                                             

% Grab your pick of algorithms.
%algorithms = {twf,plamp,plift};
algorithms = {wf,twf,taf,rwf,raf,fienup,gs,ampflow,kac,cd,scgm,pmax,plamp,plift};
algorithms = {pmax};

%% 2. Run benchmark
benchmarkTransmissionMatrix(imageSize, datasetSelection, residualConstant, algorithms)
