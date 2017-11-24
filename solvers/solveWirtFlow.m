%                           solveWirtFlow.m
%
%  Implementation of the Wirtinger Flow (WF) algorithm proposed in the
%  paper.
%
%% I/O
%  Inputs:
%     A:    m x n matrix or a function handle to a method that
%           returns A*x.     
%     At:   The adjoint (transpose) of 'A'. If 'A' is a function handle,
%           'At' must be provided.
%     b0:   m x 1 real,non-negative vector consists of all the measurements.
%     x0:   n x 1 vector. It is the initial guess of the unknown signal x.
%     opts: A struct consists of the options for the algorithm. For details,
%           see header in solvePhaseRetrieval.m or the User Guide.
%
%     Note: When a function handle is used, the
%     value of 'At' (a function handle for the adjoint of 'A') must be 
%     supplied.
% 
%  Outputs:
%     sol:  n x 1 vector. It is the estimated signal.
%     outs: A struct consists of the convergence info. For details,
%           see header in solvePhaseRetrieval.m or the User Guide.
%  
%  
%  See the script 'testWirtFlow.m' for an example of proper usage of this 
%  function.
%
%% Notations
%  x is the estimation of the signal. y is the vector of measurements such
%  that yi = |<ai,x>|^2 for i = 1,...,m
%
%% Algorithm Description
%  WF successively refines the estimate via an update rule that bears a
%  strong resemblance to a gradient descent scheme. Specifically, at each
%  iteration, x = x + mu/m * gradient log-likelihood of x given y For the
%  detailed formulation of "gradient log-likelihood of x given y" and a
%  detailed explanation of the theory, see the WF paper referenced below.
%  
%% References
%  Paper Title:   Phase Retrieval via Wirtinger Flow: Theory and Algorithms
%  Place:         Chapter 2.3
%  Authors:       Emmanuel Candes, Xiaodong Li, Mahdi Soltanolkotabi
%  arXiv Address: https://arxiv.org/abs/1407.1065
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs] = solveWirtFlow(A, At, b0, x0, opts)   
    innerOpts = struct;
    innerOpts.maxIters = opts.maxIters;
    innerOpts.maxTime = opts.maxTime;
    innerOpts.tol = opts.tol;
    innerOpts.verbose = opts.verbose;
    innerOpts.recordTimes = opts.recordTimes;
    innerOpts.recordResiduals = opts.recordResiduals;
    innerOpts.recordMeasurementErrors = opts.recordMeasurementErrors;
    innerOpts.recordReconErrors = opts.recordReconErrors;
    innerOpts.xt = opts.xt;
    
    innerOpts.searchMethod = opts.searchMethod;
    innerOpts.betaChoice = opts.betaChoice;
    
    [sol, outs] = gradientDescentSolver(A, At, x0, b0, @updateObjective, innerOpts);
    
    function [f, gradf] = updateObjective(~, ~)
        f = @(z) 0.5 * norm(abs(z).^2 - b0.^2)^2;
        gradf = @(z) (abs(z).^2 - b0.^2) .* z;
    end
end
