%% ------------------------solveKaczmarzSimple.m--------------------------------------


% Solver for the Kaczmarz method as given in Algorithm 3 of the Kaczmarz
% paper. Refer to the userguide for a detailed usage of the package.

%  See the script 'testKaczmarzSimple.m' for an example of proper usage of
%  this function.

% PAPER TITLE:
%              Solving systems of phaseless equations via Kaczmarz methods:
%              A proof of concept study.

% ARXIV LINK:
%              https://arxiv.org/pdf/1502.01822.pdf

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
%             The kaczmarz method is an iterative method based on an
%             algebraic reconstruction technique where an intial guess is
%             chosen and the next iterate is obtained by projecting the
%             current iterate onto the hyperplane <a_r, x> = y_r. a_r is
%             the rth measurement row of the matrix A.

% METHOD:
%         1.) The method is described in algorithm 3 in the paper.
 

% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs] = solveKaczmarzSimple(A, At, b0, x0, opts)
    
    validateOptions(opts);    % Check the validity of algorithm-specific options
    
    % Initialization
    m = length(b0); % num of measurements
    sol = x0;
    maxDiff = -inf;

    % Initialize values potentially computed at each round.
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];
    
    
    
    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);
    
    
    startTime = tic; % Start timer
    for iter = 1 : opts.maxIters
        switch lower(opts.indexChoice)
            case 'random'
                index = randi(m);
            case 'cyclic'
                index = mod(iter-1, m) + 1;
        end
        
        % Obtain measurement vector a_i as column vector
        a = At(double(1:m == index)');
        y = b0(index);
        product = dot(a, sol);
        
        % Calculate new estimate
        newSol = sol + ((y * sign(product) - product) / norm(a)^2) * a;
        
        diff = norm(newSol - sol);
        sol = newSol;
        maxDiff = max(diff, maxDiff);

        
        % Record convergence information and check stopping condition
        % If xt is provided, reconstruction error will be computed and used for stopping
        % condition. Otherwise, residual will be computed and used for stopping
        % condition.
        if ~isempty(opts.xt)
            x = sol;
            xt = opts.xt;
            %  Compute optimal rotation
            alpha = (x(:)'*xt(:))/(x(:)'*x(:));
            x = alpha*x;
            currentReconError = norm(x-xt)/norm(xt);
            if opts.recordReconErrors
                reconErrors(iter) = currentReconError;
            end
        end

        if isempty(opts.xt) | opts.recordResiduals
            currentResid = diff / max(maxDiff, 1.0e-30);
        end

        if opts.recordResiduals
            residuals(iter) = currentResid;
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

        %  Test stopping criteria. 
        if stopNow(opts, currentTime, currentResid, currentReconError)
            break;
        end
    end
    
    % Create output according to the options chosen by user
    outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals);

    % Display verbose output if specified
    if opts.verbose == 1
        displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
    end
end

% Validate algorithm-specific options
function validateOptions(opts)
    validIndexChoices = {'cyclic', 'random'};
    checkIfInList('indexChoice',opts.indexChoice,validIndexChoices);
end
