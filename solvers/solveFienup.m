%                           solveFienup.m
%
%  Solver for Fienup algorithm.
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
%  See the script 'testFienup.m' for an example of proper usage of this 
%  function.
%
%% Notations
%  The notations mainly follow those used in Section 2 of the Fienup paper.
%  gk:    g_k   the guess to the signal before the k th round
%  gkp:   g_k'  the approximation to the signal after the k th round of 
%         iteration
%  gknew: g_k+1 the guess to the signal before the k+1 th round
%  Gkp:   G_k'  the approximation to fourier transfor of the signal after 
%               satisfying constraints on fourier-domain
%  beta:  \beta the Tuning parameter for object-domain update
%  
%% Algorithm Description
%  Fienup Algorithm is the same as Gerchberg-Saxton Algorithm except when
%  the signal is real and non-negative (or has constraint in general). When
%  this happens, the update on the object domain is different.
%  
%  Like Gerchberg-Saxton, Fienup transforms back and forth between the two
%  domains, satisfying the constraints in one before returning to the
%  other. The method has four steps (1) Left multipy the current estimation
%  x by the measurement matrix A and get Ax. (2) Keep phase, update the
%  magnitude using the measurements b0, z = b0.*sign(Ax). (3) Solve the
%  least-squares problem
%           sol = \argmin ||Ax-z||^2
%      to get our new estimation x. We use Matlab built-in solver lsqr()
%      for this least square problem.
%  (4) Impose temporal constraints on x(This step is ignored when there is
%  no constraints)
%
%  For a detailed explanation, see the Fienup paper referenced below.
%

%% References
%  Paper Title:   Phase retrieval algorithms: a comparison
%  Place:         Section II for notation and Section V for the
%                 Input-Output Algorithm
%  Authors:       J. R. Fienup
%  Address: https://www.osapublishing.org/ao/abstract.cfm?uri=ao-21-15-2758
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 

function [sol,outs] = solveFienup(A,At,b0,x0,opts)

    validateOptions(opts);    % Check the validity of algorithm-specific options

    % Initialization
    gk = x0;                      % Initial Guess, corresponds to g_k in the paper
    gkp = x0;                     % corresponds to g_k' in the paper
    gknew = x0;                   % corresponds to g_k+1 in the paper
    beta = opts.FienupTuning;     % GS tuning parameter

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
    
    for iter = 1 : opts.maxIters
        
        Ax = A(gk);            % Intermediate value to save repetitive computation
        Gkp = b0.*sign(Ax);    % Calculate the initial spectral magnitude, G_k' in the paper.    



        % -----------------------------------------------------------------------
        % Record convergence information and check stopping condition
        % If xt is provided, reconstruction error will be computed and used for stopping
        % condition. Otherwise, residual will be computed and used for stopping
        % condition.
        if ~isempty(opts.xt)
            x = gk;
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
            currentResid = norm(At(Ax-Gkp))/norm(Gkp);
        end

        if opts.recordResiduals
            residuals(iter) = currentResid;
        end

        currentTime = toc(startTime);                % Record elapsed time so far
        if opts.recordTimes
            solveTimes(iter) = currentTime;
        end
        if opts.recordMeasurementErrors
            currentMeasurementError = norm(abs(A(gk)) - b0) / norm(b0);
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
        % -----------------------------------------------------------------------



        % Solve the least-squares problem 
        %  gkp = \argmin ||Ax-Gkp||^2.
        % If A is a matrix,
        %  gkp = inv(A)*Gkp
        % If A is a fourier transform( and measurements are not oversampled i.e. m==n),
        %  gkp = inverse fourier transform of Gkp  
        % Use the evalc() to capture text output, thus preventing
        % the conjugate gradient solver from printing to the screen.
        evalc('gkp=lsqr(@Afun,Gkp,opts.tol/100,opts.maxInnerIters,[],[],gk)');

        % If the signal is real and non-negative, Fienup updates object domain
        % following the constraint
        if opts.isComplex == false & opts.isNonNegativeOnly == true
             
            inds = gkp<0;  % Get indices that are outside the non-negative constraints
                           % May also need to check if isreal
            inds2 = ~inds; % Get the complementary indices

            % hybrid input-output (see Section V, Equation (44))
            gknew(inds) = gk(inds) - beta*gkp(inds);
            gknew(inds2) = gkp(inds2);
        else % Otherwise, its update is the same as the GerchBerg-Saxton algorithm
            gknew = gkp;
        end
        gk = gknew;               % update gk
        
    end

    sol = gk;
    % Create output according to the options chosen by user
    outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals);

    % Display verbose output if specified
    if opts.verbose == 1
        displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
    end
end


% Check the validify of algorithm specific options
function validateOptions(opts)
    checkIfNumber('FienupTuning', opts.FienupTuning);
end
