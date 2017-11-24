%%                   testGerchbergSaxton.m

% This test file implements the Gerchberg-Saxton solver. The code builds
% a synthetic formulation of the Phase Retrieval problem, b = |Ax| and
% computes an estimate to x. The code finally plots a convergence curve
% and also makes a scatter plot of the true vs recovered solution.

% PAPER TITLE:
%              A practical algorithm for the determination of the phase
%              from image and diffraction plane pictures

% ARXIV LINK:
%              http://www.u.arizona.edu/~ppoon/GerchbergandSaxton1972.pdf


% 1) Each test script starts out by defining the length of the unknown
% signal, n and the number of measurements, m. These mesurements can be
% made complex by setting the isComplex flag to be true.
%
% 2) We then build the test problem by invoking the function
% 'buildTestProblem' which generates random gaussian measurements according
% to the user's choices in step(1). The function returns the measurement 
% matrix 'A', the true signal 'xt' and the measurements 'b0'.
%
% 3) We set the options for the PR solver. For example, the maximum
% number of iterations, the tolerance value, the algorithm and initializer
% of choice. These options are controlled by setting the corresponding
% entries in the 'opts' struct.  Please see the user guide for a complete 
% list of options.
%
% 4) We solve the phase retrieval problem by running the following line 
% of code:
%   >>  [x, outs, opts] = solvePhaseRetrieval(A, A', b0, n, opts)
% This solves the problem using the algorithm and initialization scheme
% specified by the user in the struct 'opts'.
%
% 5) Determine the optimal phase rotation so that the recovered solution
% matches the true solution as well as possible.
%
% 6) Report the relative reconstruction error. Plot residuals (a measure
% of error) against the number of iterations and plot the real part of the 
% recovered signal against the real part of the original signal.


% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017
%% -----------------------------START----------------------------------


clc
clear
close all

n = 100;          % Dimension of unknown vector
m = 5 * n;        % Number of measurements
isComplex = true; % If the signal and measurements are complex

%%  Build a random test problem
fprintf('Building test problem...\n');
[A, xt, b0] = buildTestProblem(m, n, isComplex);

% Options
opts = struct;
opts.initMethod = 'Truncatedspectral';
opts.algorithm = 'GerchbergSaxton';
opts.isComplex = isComplex;
opts.maxInnerIters = 500;
opts.tol = 1e-6;
opts.verbose = 2;

%% Try to recover x
fprintf('Running algorithm...\n');
[x, outs, opts] = solvePhaseRetrieval(A, A', b0, n, opts);

%% Determine the optimal phase rotation so that the recovered solution
%  matches the true solution as well as possible.  
alpha = (x'*xt)/(x'*x);
x = alpha * x;

%% Determine the relative reconstruction error.  If the true signal was 
%  recovered, the error should be very small - on the order of the numerical
%  accuracy of the solver.
reconError = norm(xt-x)/norm(xt);
fprintf('relative recon error = %d\n', reconError);

% Plot a graph of error(definition depends on if opts.xt is provided) versus
% the number of iterations.
plotErrorConvergence(outs, opts)

% Plot a graph of the recovered signal x against the true signal xt.
plotRecoveredVSOriginal(x,xt);
