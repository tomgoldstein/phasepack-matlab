%% ------------------------solveRWF.m--------------------------------------


% Solver for Reweighted Wirtinger Flow as given in Algorithm 1
% of the Reweightde Wirtinger Flow (TAF) paper. Refer to the userguide for
% a detailed usage of the package.

% PAPER TITLE:
%              Phase Retrieval via Reweighted Wirtinger Flow
% ARXIV LINK:
%              https://doi.org/10.1364/AO.56.002418

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
%             a weighted truncated iterative descent method that removes
%             those measurement vector a_m that are not correlated to the
%             intial guess. This truncation and weighting step is done at
%             each iteration of the gradient descent method.

% METHOD:
%         1) Our implementation uses FASTA, a fast gradient solver.
%
%         2) Set the objectve f = @(z) 1/2 * sum(weights .* 
%             (abs(z) - b0).^2). The weights are computed according to
%             equation (14) in algorithm 1 of the proposing paper
%
%         3) The gradient is grad = @(z) weights .* (z - b0 .* sign(z)).
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs] = solveRWF(A, At, b0, x0, opts)
%     addpath('solvers/linesearch');
    
    m = length(b0);
    
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
    
    innerOpts.updateObjectivePeriod = opts.reweightPeriod;
    innerOpts.searchMethod = opts.searchMethod;
    innerOpts.betaChoice = opts.betaChoice;
    
    [sol, outs] = gradientDescentSolver(A, At, x0, b0, @updateObjective, innerOpts);
    
    function [f, gradf] = updateObjective(~, Ax)
        weights = 1 ./ (abs(abs(Ax).^2 - b0.^2) + opts.eta*ones(m, 1));
        s = sum(weights);
        f = @(z) 0.5/s * sum(weights .* (abs(z).^2 - b0.^2).^2);
        gradf = @(z) (1.0/s)*weights .* (abs(z).^2 - b0.^2) .* z;
    end
end