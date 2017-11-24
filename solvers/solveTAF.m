%% ------------------------solveTAF.m--------------------------------------


% Solver for Truncated Amplitude Flow as given in Algorithm 1
% of the Truncated Amplitude Flow (TAF) paper.Refer to the userguide for
% a detailed usage of the package.

%  See the script 'testTAF.m' for an example of proper usage of
%  this function.

% PAPER TITLE:
%              Solving Systems of Random Quadratic Equations via Truncated
%              Amplitude Flow.

% ARXIV LINK:
%              https://arxiv.org/pdf/1605.08285.pdf

% INPUTS:
%         A:   Function handle/numerical matrix for data matrix A. The rows
%              of this matrix are the measurement vectors that produce
%              amplitude measurements '\psi'.
%         At:  Function handle/numerical matrix for A transpose.
%         b0:  Observed data vector consisting of amplitude measurements
%              generated from b0 = |A*x|. We assign it to 'psi' to be
%              consistent with the notation in the paper.
%         x0:  The initial vector to be used by any solver. 
%        opts: struct consists of the options for the algorithm. For
%              details,see header in solvePhaseRetrieval.m or the User
%              Guide.

% OUPTUT :
%         sol: n x 1 vector. It is the estimated signal.
%        outs: A struct consists of the convergence info. For details,
%              see header in solvePhaseRetrieval.m or the User Guide.

% Note:        When a function handle is used, the value of 'n' (the length
%              of the unknown signal) and 'At' (a function handle for the
%              adjoint of 'A') must be supplied. When 'A' is numeric, the
%              values of 'At' and 'n' are ignored and inferred from the
%              arguments


% DESCRIPTION:
%             The spectral method proposed by Candes et al. in their
%             seminal Wirtinger Flow paper suffers from fat tails that can
%             distort the true solution. A possible way to overcome this is
%             to truncate the outliers. In this paper, the authors propose
%             a truncated iterative descent method that removes those
%             measurement vector a_m that are not orthogonal to the intial
%             guess. This truncationn step is done at each iteration of the
%             gradient descent method.

% METHOD:
%         1) Our implementation uses FASTA, a fast gradient solver.
%            
%         2) Set the objectve f = 1/2 * sum(mask .* (abs(z) - b0).^2).
%             The mask contains the truncated vectors that need to be
%             included.
%
%         3) The gradient is grad = mask .* (z - b0 .* sign(z)). 
%     
%         4) The truncation is only updated periodically to enable the
%         solver to use fast adaptive stepsize rules.

% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs] = solveTAF(A, At, b0, x0, opts)
%     addpath('solvers/linesearch');
    
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
    
    innerOpts.updateObjectivePeriod = opts.truncationPeriod;
    innerOpts.searchMethod = opts.searchMethod;
    innerOpts.betaChoice = opts.betaChoice;
    
    [sol, outs] = gradientDescentSolver(A, At, x0, b0, @updateObjective, innerOpts);
    
    function [f, gradf] = updateObjective(~, Ax)
        mask = abs(Ax) >= 1/(1+opts.gamma) * b0;
        s = sum(mask);
        f = @(z) (0.5/s) * sum(mask .* (abs(z) - b0).^2);
        gradf = @(z) (1.0/s) * mask .* (z - b0 .* sign(z));
    end
end