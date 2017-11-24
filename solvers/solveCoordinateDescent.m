%                           solveCoordinateDescent.m
%
%  Implementation of the Coordinate Descent algorithm proposed in the paper.
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
%  See the script 'testCoordinateDescent.m' for an example of proper usage
%  of this function.
%
%% Notations
%  Notations mainly follow the paper's notation.
%
%% Algorithm Description
%  CD is an iterative procedure that successively minimizes the objective
%  function along coordinate directions. A single unknown is solved at each
%  iteration while all other variables are kept fixed. As a result, only
%  minimization of a univariate quartic polynomial is needed which is
%  easily achieved by finding the closed-form roots of a cubic equation.
%  
%  Specifically, the method has the following steps: It keeps running until
%  the normalized gradient becomes smaller than the tolerance (1) At each
%  iteration, use the selected rule to choose an index i. (2) Minimize the
%  objective f with respect to the ith variable while keeping
%      all other 2n-1 (both real and imaginary parts are variables so there
%      are 2n in total) variables fixed by solving the cubic equation to
%      get alpha that minimize the objective along the ith variable.
%  (3) update the estimate along the ith variable by alpha.
%   
%% References
%  Paper Title:   Coordinate Descent Algorithms for Phase Retrieval
%  Place:         Chapter II.B
%  Authors:       Wen-Jun Zeng, H. C. So
%  arXiv Address: https://arxiv.org/abs/1706.03474
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 

function [sol, outs] = solveCoordinateDescent(A, At, b0, x0, opts)
    
    validateOptions(opts); % Check the validity of algorithm-specific options
    
    % Initialization
    m = length(b0);
    n = length(x0);
    sol = x0;
    Ax = A(sol);

    % Initialize values potentially computed at each round.
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];
    
    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);

    maxDiff = -inf;
    
    C = zeros(m, 3);
    
    % Used to compute gradient of objective function (used by greedy index
    % choice rule)
    f = @(z) (abs(z).^2 - b0.^2) .* z;
    
    
    startTime = tic; % Begin timer
    for iter = 1 : opts.maxIters 
        switch lower(opts.indexChoice)
            case 'random'
                if opts.isComplex
                    index = randi(2*n);
                else
                    index = randi(n);
                end
            case 'cyclic'
                if opts.isComplex
                    index = mod(iter-1, 2*n) + 1;
                else
                    index = mod(iter-1, n) + 1;
                end
            case 'greedy'
                grad = At(f(A(sol)));
                
                if opts.isComplex
                    grad_bar = [real(grad); imag(grad)];
                else
                    grad_bar = grad;
                end
                [~, index] = max(abs(grad_bar));
        end
        
        vals = conj(A(double(1:n == mod(index-1, n)+1)'));
        for j = 1 : m 
            if index > n
                C(j, 2) = 2*imag(Ax(j) * vals(j));
            else
                C(j, 2) = 2*real(Ax(j) * vals(j));
            end
            C(j, 3) = abs(vals(j))^2;
            C(j, 1) = abs(Ax(j))^2;
        end
        
        d_4 = sum(C(:, 3).^2);
        d_3 = sum(2 * C(:, 3) .* C(:, 2));
        d_2 = sum(C(:, 2).^2 + 2 * C(:, 3) .* (C(:, 1) - b0.^2));
        d_1 = sum(2 * C(:, 2) .* (C(:, 1) - b0.^2));
        
        % Desired alpha is a root of polynomial
        alphas = roots([4*d_4 3*d_3 2*d_2 d_1])';
        % Select only the real roots
        alphas = alphas(imag(alphas) == 0);
        % Function of alpha to be minimized
        g = @(x) d_4*x.^4 + d_3*x.^3 + d_2*x.^2 + d_1*x;
        % Find index of best alpha
        [~, idx] = min(g(alphas));
        alpha = alphas(idx);
        
        % Update x
        if (index > n)
            a_j = A(double(1:n == index - n)');
            sol(index - n) = sol(index - n) + 1i * alpha;
            Ax = Ax + (1i * alpha) * a_j;
        else
            a_j = A(double(1:n == index)');
            sol(index) = sol(index) + alpha;
            Ax = Ax + alpha * a_j;
        end
        
        diff = abs(alpha);
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
    validIndexChoices = {'cyclic', 'random', 'greedy'}; 
    checkIfInList('indexChoice',opts.indexChoice,validIndexChoices);
end
