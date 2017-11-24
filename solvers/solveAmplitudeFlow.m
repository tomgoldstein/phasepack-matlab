%% ------------------------solveAmplitudeFlow.m----------------------------


% Solver for Amplitude Flow as given in Algorithm 1 of the Truncated
% Amplitude Flow (TAF) paper minus the truncations.  The idea was to test the
% amplitude based objective function without removing any outliers. The
% advantage is that we can easily call FASTA. Refer to the userguide for a
% detailed usage of the package.

%  See the script 'testAmplitudeFlow.m' for an example of proper usage of
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
%             The wirtinger flow algorithm uses the squared magnitude
%             objective. This generates a slightly different gradient then
%             the amplitude based objective which is not squared. The TAF
%             paper proposes an amplitude based objective optimization
%             using truncation. This code implements a plain vanilla
%             gradient descent without truncation.

% METHOD:
%         1.) Our implementation uses FASTA, A Forward Backward Splitting
%             Package. FASTA can solve gradient descent schemes provided
%             the objective function expression, its gradient expression,
%             the non smooth term and its proximal operator
%

%         2.) Set the objectve f = @(z) 1/2 * norm(abs(z) - b0)^2.
%             The mask contains the truncated vectors that need to be
%             included.

%         3.) The gradient grad = @(z) (z - b0 .* sign(z)). There is no
%             non smooth term hence g = 0. Conseqently prox = @(z) z
%     
%         4.) Send all the above parameters to FASTA which then spits out
%             the estimated signal. 

% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs] = solveAmplitudeFlow(A, At, b0, x0, opts)
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
    
    innerOpts.searchMethod = opts.searchMethod;
    innerOpts.betaChoice = opts.betaChoice;
    
    [sol, outs] = gradientDescentSolver(A, At, x0, b0, @updateObjective, innerOpts);
    
    function [f, gradf] = updateObjective(~, ~)
        f = @(z) 1/2 * norm(abs(z) - b0)^2;
        gradf = @(z) (z - b0 .* sign(z));
    end
end