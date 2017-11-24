%                           solvePhaseMax.m
%
%  Implementation of the PhaseMax algorithm proposed in the paper using
%  FASTA. Note: For this code to run, the solver "fasta.m" must be in your
%  path.
%
%% I/O
%  Inputs:
%     A:    m x n matrix or a function handle to a method that
%           returns A*x.     
%     At:   The adjoint (transpose) of 'A'. If 'A' is a function handle, 'At'
%           must be provided.
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
%  See the script 'testPhaseMaxGaussian.m' and 'testPhaseMaxFourier.m' for 
%  two examples of proper usage of this function.
%
%% Notations
%  
%  
%% Algorithm Description
%  Solve the PhaseMax signal reconstruction problem
%         maximize <x0,x>
%         subject to |Ax|<=b0
% 
%  The problem is solved by approximately enforcing the constraints using a
%  quadratic barrier function.  A continuation method is used to increase
%  the strength of the barrier function until a high level of numerical
%  accuracy is reached.
%  The objective plus the quadratic barrier has the form
%     <-x0,x> + 0.5*max{|Ax|-b,0}^2.
%
%  For a detailed explanation, see the PhaseMax paper referenced below. For
%  more details about FASTA, see the FASTA user guide, or the paper "A
%  field guide to forward-backward splitting with a FASTA implementation."
%
%% References
%  Paper Title:   PhaseMax: Convex Phase Retrieval via Basis Pursuit
%  Authors:       Tom Goldstein, Christoph Studer
%  arXiv Address: https://arxiv.org/abs/1610.07531
%  
%  Copyright Goldstein & Studer, 2016.  For more details, visit 
%  https://www.cs.umd.edu/~tomg/projects/phasemax/

function [sol, outs] = solvePhaseMax(A,At,b0,x0,opts)
    % Initialization
    m = length(b0);              % number of measurements
    n = length(x0);              % length of the unknown signal
    remainIters = opts.maxIters; % The remaining fasta iterations we have. 
                                 % It's initialized to opts.maxIters.

    %  Normalize the initial guess relative to the number of measurements
    x0 = (x0/norm(x0(:)))*mean(b0(:))*(m/n)*100;

    %  re-scale the initial guess so that it approximately satisfies |Ax|=b
    sol = x0.* min(b0./abs(A(x0)));

    % Initialize values potentially computed at each round.
    ending = 0;       % Indicate whether any of the ending condition (except maxIters) has been met in FASTA.
    iter = 0;
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];

    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);

    %%  Define objective function components for the gradient descent method FASTA
    f = @(z) 0.5*norm(max(abs(z)-b0,0))^2;     % f(z) = 0.5*max{|x|-b,0}^2 : This is the quadratic penalty function
    gradf = @(z)  (sign(z).*max(abs(z)-b0,0)); % The gradient of the quadratic penalty
    
    % Options to hand to fasta
    fastaOpts.maxIters = opts.maxIters;  
    fastaOpts.stopNow = @(x, iter, resid, normResid, maxResid, opts) ...
        processIteration(x, resid);  % Use customized stopNow in order to get
                                         % solveTime, residual and error at each iteration.
    fastaOpts.verbose=0;
    startTime = tic;                    % Start timer
    constraintError = norm(abs(A(sol))-b0);  % Keep track of the current error in the solution.
    while remainIters > 0 & ~ending     % Iterate over continuation steps
        g = @(x) -real(x0'*x);          % The linear part of the objective
        proxg = @(x,t) x+t*x0;          % The proximal operator of the linear objective
        fastaOpts.tol = norm(x0)/100;   % use a tighter tolerance when the solution is more exact
        % Call FASTA to solve the inner minimization problem
        [sol, fastaOuts] = fasta(A, At, f, gradf, g, proxg, sol, fastaOpts);  % Call a solver to minimize the quadratic barrier problem

        fastaOpts.tau = fastaOuts.stepsizes(end);     % Record the most recent stepsize for recycling.
        x0 = x0/10;                                   % do continuation - this makes the quadratic penalty stronger

        % Update the max number of iterations for fasta
        remainIters = remainIters - fastaOuts.iterationCount;
        fastaOpts.maxIters = min(opts.maxIters, remainIters);
        
        % Monitor convergence and check stopping conditions
        newConstraintError = norm(max(abs(A(sol))-b0,0));
        relativeChange = abs(constraintError-newConstraintError)/norm(b0);
        if relativeChange <opts.tol   % terminate when error reduction stalls
            break;
        end
        constraintError = newConstraintError;
    end


    % Create output according to the options chosen by user
    outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals);

    % Display verbose output if specified
    if opts.verbose == 1
        displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
    end


    % Runs code upon each FASTA iteration. Returns whether FASTA should
    % terminate.
    function stop = processIteration(x, residual)
        iter = iter + 1;
        % Record convergence information and check stopping condition
        % If xt is provided, reconstruction error will be computed and used for stopping
        % condition. Otherwise, residual will be computed and used for stopping
        % condition.
        if ~isempty(opts.xt)
            xt = opts.xt;
            %  Compute optimal rotation
            alpha = (x(:)'*xt(:))/(x(:)'*x(:));
            x = alpha*x;
            currentReconError = norm(x-xt)/norm(xt);
            if opts.recordReconErrors
                reconErrors(iter) = currentReconError;
            end
        end

        if isempty(opts.xt)
            currentResid = residual;
        end

        if opts.recordResiduals
            residuals(iter) = residual;
        end

        currentTime = toc(startTime);                % Record elapsed time so far
        if opts.recordTimes
            solveTimes(iter) = currentTime;
        end
        if opts.recordMeasurementErrors
            currentMeasurementError = norm(abs(A(sol)) - b0) / norm(b0);
            measurementErrors(iter) = currentMeasurementError;
        end
       
        % Display verbose output if specified
        if opts.verbose == 2
            displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
        end

        % Test stopping criteria.
        stop = false;
        if currentTime >= opts.maxTime % Stop if we're run over the max runtime
            stop = true;
        end
        if ~isempty(opts.xt) % If true solution is specified, terminate when close to true solution
            assert(~isempty(currentReconError),'If xt is provided, currentReconError must be provided.');
            stop = stop || currentReconError < opts.tol;
            ending = stop;  % When true, this flag will terminate outer loop
        end
        stop = stop || residual < fastaOpts.tol;    % Stop FASTA is the tolerance is reached
    end

end

