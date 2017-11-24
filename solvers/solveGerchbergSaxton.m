%                           solveGerchbergSaxton.m
%
%  Solver for Gerchberg-Saxton algorithm.
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
%  See the script 'testGerchbergSaxton.m' for an example of proper usage of
%  this function.
%
%% Notations
%  z is the spectral magnitude on the fourier domain(z = b0.*sign(Ax)). 
%  It has the phase of the fourier transform Ax of the estimation
%  and the magnitude of the measurements b0.
%  
%% Algorithm Description
%  One simply transforms back and forth between the two domains, satisfying
%  the constraints in one before returning to the other. The method has
%  three steps 
% (1) Left multipy the current estimation x by the measurement
%  matrix A and get Ax. 
% (2) Keep phase, update the magnitude using the
%  measurements b0, z = b0.*sign(Ax). 
% (3) Solve the least-squares problem
%           sol = \argmin ||Ax-z||^2
%      to get our new estimation x. We use Matlab built-in solver lsqr()
%      for this least square problem.
%  (4) Impose temporal constraints on x(This step is ignored now since we
%      don't assume constraints on our unknown siganl)
%
%  For a detailed explanation, see the paper referenced below.
%

%% References
%  The paper that this implementation follows
%  Paper Title:   Fourier Phase Retrieval: Uniqueness and Algorithms
%  Place:         Chapter 4.1 Alternating projection algorithms, Algorithm 1
%  Authors:       Tamir Bendory, Robert Beinert, Yonina C. Eldar
%  arXiv Address: https://arxiv.org/abs/1705.09590
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 

function [sol,outs] = solveGerchbergSaxton(A,At,b0,x0,opts)

    validateOptions(opts);    % Check the validity of algorithm-specific options

    % Initialization
    sol = x0;                 % Initial Guess

    % Initialize values potentially computed at each round.
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];
 
    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);

    % Build a function handle for matlab's conjugate-gradient solver
    function y = Afun(x,transp_flag)
       if strcmp(transp_flag,'transp')       % y = A'*x
          y = At(x);
       elseif strcmp(transp_flag,'notransp') % y = A*x
          y = A(x);
       end
    end
    
    startTime = tic; % Start timer
    z = b0.*sign(A(sol)); % Calculate the initial spectral magnitude.
    for iter = 1 : opts.maxIters
        % Solve the least-squares problem 
        %  sol = \argmin ||Ax-z||^2.
        % If A is a matrix,
        %  sol = inv(A)*z
        % If A is a fourier transform( and measurements are not oversampled i.e. m==n),
        %  sol = inverse fourier transform of z  
        % Use the evalc() to capture text output, thus preventing
        % the conjugate gradient solver from printing to the screen.
        evalc('sol=lsqr(@Afun,z,opts.tol/100,opts.maxInnerIters,[],[],sol)');
        Ax = A(sol);              % Intermediate value to save repetitive computation
        z = b0.*sign(Ax);         % Update the spectral magnitude

       
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
        
        if isempty(opts.xt) || opts.recordResiduals
            currentResid = norm(At(Ax-z))/norm(z);
        end
       
        if opts.recordResiduals
            residuals(iter) = currentResid;
        end

        currentTime = toc(startTime);                % Record elapsed time so far
        if opts.recordTimes
            solveTimes(iter) = currentTime;
        end
        if opts.recordMeasurementErrors
            currentMeasurementError = norm(abs(Ax) - b0) / norm(b0);
            measurementErrors(iter) = currentMeasurementError;
        end
       
        % Display verbose output if specified
        if opts.verbose == 2
            displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
        end

        % Test stopping criteria. 
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

% Check the validify of algorithm specific options
function validateOptions(opts)
    checkIfNumber('maxInnerIters', opts.maxInnerIters);
end
